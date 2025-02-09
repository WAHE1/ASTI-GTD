import dendropy
import linkedMSC
import argparse
import numpy as np
import random

# replace elements
parser = argparse.ArgumentParser()
parser.add_argument('-t','--tree', default = 'species_tree.tre')
parser.add_argument('-r','--rate', type = float, default = 1e6)
parser.add_argument('-n','--num', type = int, default = 200)
parser.add_argument('-o','--outp', default = 'gene_trees.tre')
parser.add_argument('--seed', type=int, help="Random seed for reproducibility.", default=None)

args = parser.parse_args()
tree = args.tree
rate = args.rate
num = args.num
outp = args.outp
seed = args.seed

# Set the seed if provided
if seed is not None:
    print(f"Random seed set to: {seed}")
else:
    print("No seed provided. Using default random behavior.")
rng = random.Random()
rng.seed(seed)

# read the species tree
species_tree = dendropy.Tree.get(path = str(tree),schema = "newick")

# --------------------------------------------------------------------------------------------
# generate gene to species map
num_contained = 1
gene_to_species_map = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
        containing_taxon_namespace=species_tree.taxon_namespace,
        num_contained=num_contained)

# --------------------------------------------------------------------------------------------
# write out the gene to species map for ASTRAL use
f = open("gene_to_species_map.txt","w")

species = list(species_tree.taxon_namespace)
individuals = list(gene_to_species_map)

num_species = len(species)
individuals_start_index = 0

if type(num_contained) == int:
   num_contained = [num_contained] * num_species

for i in range(num_species):
   f.write(species[i].label + ':')
   num_individuals = num_contained[i]
   if num_individuals == 1:
       f.write(individuals[individuals_start_index].label.replace(" ", "_"))
   else:
       for j in range(num_individuals - 1):
           f.write(individuals[individuals_start_index + j].label.replace(" ", "_") + ',')
       f.write(individuals[individuals_start_index + j + 1].label.replace(" ", "_"))
   individuals_start_index = individuals_start_index + num_individuals
   f.write('\n')

f.close()

# --------------------------------------------------------------------------------------------
# generate gene trees and store in a file
from dendropy import TreeList
gene_trees = TreeList()
guide_tree= dendropy.model.coalescent.contained_coalescent_tree(containing_tree=species_tree,
                                                                gene_to_containing_taxon_map=gene_to_species_map,
                                                                rng=rng)
gene_trees.append(guide_tree)
for i in range(num-1):
    # Set different seed for each gene tree generation to avoid a set of same trees
    gene_tree = linkedMSC.multispecies_linked_coalescent(species_tree, guide_tree, recombination_rate=rate, pop_size=1, rng=rng)
    gene_trees.append(gene_tree)
    guide_tree = gene_tree

## To make taxon lable from taxon_1 to taxon
#for i in range(num):
#    for t in gene_trees[i].leaf_nodes():
#        t.taxon.label = t.taxon.label.split( )[0]

gene_trees.write(path=outp,schema='newick')
