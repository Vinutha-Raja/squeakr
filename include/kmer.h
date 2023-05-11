/*
 * ============================================================================
 *
 *        Authors:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 *                  Rob Johnson <robj@vmware.com>   
 *                  Rob Patro (rob.patro@cs.stonybrook.edu)
 *
 * ============================================================================
 */

#ifndef _KMER_H_
#define _KMER_H_

#include <stdio.h>
#include <string>
#include <unordered_set>

#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) \
												- 1ULL)
enum DNA_MAP {C, A, T, G};  // A=1, C=0, T=2, G=3

class Kmer {
	public:
		static char map_int(uint8_t base);
		static uint8_t map_base(char base);
		static __int128_t str_to_int(std::string str);
		static std::string int_to_str(__int128_t kmer, uint64_t kmer_size);
        static std::string int_to_str(__int128_t kmer, uint64_t kmer_size, std::string &res);

		static std::string int_to_str_lmer(__int128_t kmer, uint64_t kmer_size);
		static int reverse_complement_base(int x);
		static __int128_t reverse_complement(__int128_t kmer, uint64_t kmer_size);
		static bool compare_kmers(__int128_t kmer, __int128_t kmer_rev);

		static void parse_kmers(const char *filename, uint64_t ksize,
														std::unordered_set<uint64_t>& kmerset);

	private:
		Kmer();
};

#endif
