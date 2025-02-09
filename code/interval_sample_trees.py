import dendropy
import numpy as np
import argparse
import itertools

# replace elements
parser = argparse.ArgumentParser()
parser.add_argument('-i','--intp', default = 'gene_trees.tre')
parser.add_argument('-s','--step', type = int, default = 1)
parser.add_argument('-o','--outp', default = 'sampledTrees.tre')

args = parser.parse_args()
intp = args.intp
step = args.step

trees = dendropy.TreeList.get(path = str(intp),schema = "newick")
n = len(trees)

indices = []
flag = 0
while flag < n:
    indices.append(flag)
    flag = flag + step
len(indices)

sampledTrees = dendropy.TreeList()
for i in indices:
    sampledTrees.append(trees[i])

sampledTrees.write(path = outp, schema = "newick")
