/*
 * ============================================================================
 *
 *        Authors:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 *                  Rob Johnson <robj@vmware.com>   
 *                  Rob Patro (rob.patro@cs.stonybrook.edu)
 *
 * ============================================================================
 */

#include <iostream>
#include <algorithm>
#include <cstring>
#include <vector>
#include <set>
#include <unordered_set>
#include <bitset>
#include <cassert>
#include <fstream>

#include <boost/thread/thread.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include <boost/atomic.hpp>

#include <time.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>

#include "clipp.h"
#include "ProgOpts.h"
#include "SqueakrFS.h"
#include "squeakrconfig.h"
#include "gqf_cpp.h"
#include "chunk.h"
#include "kmer.h"
#include "reader.h"
#include "util.h"
#include "minimizer.h"

#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) \
												- 1ULL)
#define QBITS_LOCAL_QF 16
#define MAX_NUM_THREADS 64

typedef struct {
	CQF<KeyObject> *local_cqf;
	CQF<KeyObject> *main_cqf;

	uint32_t count {0};
	uint32_t ksize {28};
    uint32_t lsize {15};
    bool exact;
	spdlog::logger* console{nullptr};
} flush_object;

std::set<std::string> lmerSet;
/*create a multi-prod multi-cons queue for storing the chunk of fastq file.*/
boost::lockfree::queue<file_pointer*, boost::lockfree::fixed_sized<true> > ip_files(64);
boost::atomic<int> num_files {0};

/* dump the contents of a local QF into the main QF */
static bool dump_local_qf_to_main(flush_object *obj)
{
	CQF<KeyObject>::Iterator it = obj->local_cqf->begin();
	do {
		KeyObject hash = it.get_cur_hash();
		int ret = obj->main_cqf->insert(hash, QF_WAIT_FOR_LOCK | QF_KEY_IS_HASH);
		if (ret == QF_NO_SPACE) {
			obj->console->error("The CQF is full. Please rerun the with a larger size.");
			return false;
		}
		++it;
	} while (!it.done());
	obj->local_cqf->reset();

	return true;
}

bool lmer_flag;
/* convert a chunk of the fastq file into kmers */
bool reads_to_kmers(chunk &c, flush_object *obj, MinimizerScanner &scanner)
{
	auto fs = c.get_reads();
	auto fe = c.get_reads();
	auto end = fs + c.get_size();
	while (fs && fs!=end) {
		fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore the first line
		fs++; // increment the pointer

		fe = static_cast<char*>(memchr(fs, '\n', end-fs)); // read the read
		std::string read(fs, fe-fs);

		if (read.length() < obj->ksize) // start with the next read if length is smaller than K
			goto next_read;
		{

            scanner.LoadSequence(read);

            uint64_t *mmp;
            while ((mmp = scanner.NextMinimizer()) != nullptr){
                lmerSet.insert(Kmer::int_to_str(*mmp, obj->lsize));
            }
		}

next_read:
		fs = ++fe;		// increment the pointer
		fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore one line
		fs++; // increment the pointer
		fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore one more line
		fs++; // increment the pointer
	}
	free(c.get_reads());

	return true;
}

/* read a part of the fastq file, parse it, convert the reads to kmers, and
 * insert them in the CQF
 */
static bool fastq_to_uint64kmers_prod(flush_object* obj)
{
	file_pointer* fp;
    MinimizerScanner scanner(obj->ksize, obj->lsize, 0, true, 0);
	while (num_files) {
		while (ip_files.pop(fp)) {
			if (reader::fastq_read_parts(fp->mode, fp)) {
				ip_files.push(fp);
				chunk c(fp->part, fp->size);
				if (!reads_to_kmers(c, obj, scanner)) {
					obj->console->error("Insertion in the CQF failed.");
					abort();
				}
			} else {
				/* close the file */
				if (fp->mode == 0)
					fclose(fp->freader->in);
				else if (fp->mode == 1)
					gzclose(fp->freader->in_gzip);
				else if (fp->mode == 2)
					if (fp->freader->in) {
						BZ2_bzReadClose(&(fp->freader->bzerror), fp->freader->in_bzip2);
						fclose(fp->freader->in);
					}
				delete[] fp->part_buffer;
				delete fp;
				num_files--;
			}
		}
	}
	if (obj->count) {
		if (!dump_local_qf_to_main(obj)) {
			obj->console->error("Insertion in the CQF failed.");
			abort();
		}
		obj->count = 0;
	}

	return true;
}


/* main method */
int count_main(CountOpts &opts)
{
	int mode = 0;

	spdlog::logger* console = opts.console.get();

	if (opts.exact && opts.ksize > 32) {
		console->error("Does not support k-mer size > 32 for squeakr-exact.");
		return 1;
	}

	if (opts.numthreads == 0)
		opts.numthreads = std::thread::hardware_concurrency();

	enum qf_hashmode hash = QF_HASH_DEFAULT;
	int num_hash_bits = opts.qbits+8;	// we use 8 bits for remainders in the main QF
	if (opts.exact) {
		num_hash_bits = 2*opts.lsize; // Each base 2 bits.
		hash = QF_HASH_INVERTIBLE;
	}

	std::string ser_ext(".squeakr");
	std::string log_ext(".log");
	std::string cluster_ext(".cluster");
	std::string freq_ext(".freq");
	struct timeval start1, start2, end1, end2;
	struct timezone tzp;
	uint32_t OVERHEAD_SIZE = 65535;

	std::string filepath(opts.filenames.front());
	auto const pos = filepath.find_last_of('.');
	std::string input_ext = filepath.substr(pos + 1);
	if (input_ext == std::string("fastq") || input_ext == std::string("fq"))
		mode = 0;
	else if (input_ext == std::string("gz"))
		mode = 1;
	else if (input_ext == std::string("bz2"))
		mode = 2;
	else {
		console->error("Does not support this input file type.");
		return 1;
	}

	for( auto& fn : opts.filenames ) {
		auto* fr = new reader;
		if (reader::getFileReader(mode, fn.c_str(), fr)) {
			file_pointer* fp = new file_pointer;
			fp->mode = mode;
			fp->freader.reset(fr);
			fp->part_buffer = new char[OVERHEAD_SIZE];
			ip_files.push(fp);
			num_files++;
		} else {
			delete fr;
		}
	}

	std::string ds_file = opts.output_file;

	boost::thread_group prod_threads;

	for (int i = 0; i < opts.numthreads; i++) {
		flush_object* obj = (flush_object*)malloc(sizeof(flush_object));
		obj->ksize = opts.ksize;
        obj->lsize = opts.lsize;
		obj->exact = opts.exact;
		obj->count = 0;
		obj->console = console;
		prod_threads.add_thread(new boost::thread(fastq_to_uint64kmers_prod,
																			obj));
	}

	console->info("Reading from the fastq file and inserting in the CQF.");
	gettimeofday(&start1, &tzp);
	prod_threads.join_all();

    ds_file.append("-set-");
    // Create an output file stream and open a file
    ds_file.append(std::to_string(opts.lsize));
    ds_file.append(".txt");
    std::ofstream outputFile(ds_file);

    // Check if the file was opened successfully
    if (!outputFile.is_open()) {
        std::cout << "Failed to open file" << std::endl;
        return 1;
    }

    for (const std::string& x : lmerSet) {
        outputFile << x << std::endl;
    }

    // Close the file
    outputFile.close();

	return 0;
}

