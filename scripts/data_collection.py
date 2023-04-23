import requests
import os

url = "https://www.cs.cmu.edu/~ckingsf/software/bloomtree/srr-list.txt"
response = requests.get(url)
text = response.text

lines = text.strip().split("\n")
SRA_list = [line.strip().split()[1] for line in lines]

# print(SRA_list)

for run in SRA_list[:50]:
    url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{}/{}/{}.fastq.gz".format(run[0:6], run, run)
    print(url)
    os.system('wget {}'.format(url))
    os.system('/home/u1306224/DS/kmer-squeakr/squeakr/squeakr count -e -k 21 -s 18 -t 1 -o /home/u1306224/DS/kmer-squeakr/squeakr/data/{}.squeakr {}.fastq.gz'.format(run, run))
    os.system('./squeakr count -e -k 21 -s 18 -l 7 -t 1 -c 10 -o data/hashsets/{} {}.fastq.gz'.format(run, run))
    os.system('./squeakr count -e -k 21 -s 18 -l 9 -t 1 -c 10 -o data/hashsets/{} {}.fastq.gz'.format(run, run))

    os.system('rm {}.fastq.gz'.format(run))
    os.system('rm wget-log*')


