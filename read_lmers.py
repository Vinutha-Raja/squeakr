import glob
from functools import cache
import seaborn as sn
import matplotlib.pyplot as plt
import os

file_sizes = {}

files = glob.glob('*-7.txt')
for file in files:
    file_sizes[file] = os.path.getsize(file)
files = sorted(file_sizes.items(), key=lambda x: x[1])
files = [f[0] for f in files]
print(files)


def create_map(file):
    doc_map = dict()
    with open(file, "r") as f:
        for line in f:
            key, value = line.strip().split(": ")
            value = int(value)
            doc_map[key] = value

    return doc_map

def create_sets(filess):
    file_sets = []
    for fi in filess:
        kmer_set = set()
        with open(fi, "r") as f:
            for line in f:
                key, value = line.strip().split(": ")
                kmer_set.add(key)
        file_sets.append(kmer_set)
    return file_sets


def weighted_jaccard(doc1, doc2):
    terms = set(doc1.keys()).intersection(set(doc2.keys()))
    min_freqs = [min(doc1.get(term, 0), doc2.get(term, 0)) for term in terms]
    numerator = sum(min_freqs)
    denominator = sum([max(doc1.get(term, 0), doc2.get(term, 0)) for term in terms])
    wjaccard = numerator / denominator
    return wjaccard

def jaccard(doc1_set, doc2_set):
    intersection_terms = doc1_set.intersection(doc2_set)
    union_terms = doc1_set.union(doc2_set)
    jc = len(intersection_terms)/len(union_terms)
    return jc

def get_maps(file_list):
    doc_list = []
    for f in file_list:
        doc_list.append(create_map(f))
    return doc_list


def find_jc(file_list, doc_list):
    data = []
    for i in range(len(file_list)):
        print("{},".format(file_list[i]), end = "")
    print("")
    for i in range(len(file_list)):
        row = []
        for j in range(0, len(file_list)):
            row.append(jaccard(doc_list[i], doc_list[j]))
            print("{},".format( jaccard(doc_list[i], doc_list[j])), end ="")
        print("")
        data.append(row)
    print(data)
    print(len(data))
    return data


docu_sets = create_sets(files)
data = find_jc(files, doc_list=docu_sets)

hm = sn.heatmap(data = data)

# displaying the plotted heatmap
plt.savefig("heat_map_7lmer.png")
plt.show()