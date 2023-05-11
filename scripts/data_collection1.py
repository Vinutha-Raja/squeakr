# Usage: python3 data_collection1.py <path of squeakr folder with master branch> <path of squeakr folder with LmerHashSet branch> <number of squeakr files> 
# Example: python3 data_collection1.py  /home/u1306224/DS/kmer-squeakr/squeakr /home/u1306224/DS/squeakr/ 2

import requests
import os
import sys

url = "https://www.cs.cmu.edu/~ckingsf/software/bloomtree/srr-list.txt"
response = requests.get(url)
text = response.text

lines = text.strip().split("\n")
SRA_list = [line.strip().split()[1] for line in lines]

print(len(SRA_list))

kmer_path = sys.argv[1]
lmer_path = sys.argv[2]
file_count = int(sys.argv[3])

for run in SRA_list[0:file_count]:
    url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{}/{}/{}.fastq.gz".format(run[0:6], run, run)
    print(url)
    os.system('wget {}'.format(url))
    os.system('{}/squeakr count -e -k 21 -s 18 -t 1 -o {}/data/{}.squeakr {}.fastq.gz'.format(kmer_path, kmer_path, run, run))
    os.system('{}/squeakr count -e -k 21 -s 18 -l 7 -t 1 -c 10 -o {}/data/{} {}.fastq.gz'.format(lmer_path, lmer_path, run, run))
    os.system('{}/squeakr count -e -k 21 -s 18 -l 9 -t 1 -c 10 -o {}/data/{} {}.fastq.gz'.format(lmer_path, lmer_path, run, run))

    os.system('rm {}.fastq.gz'.format(run))
    os.system('rm wget-log*')


