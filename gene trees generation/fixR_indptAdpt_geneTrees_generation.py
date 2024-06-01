import dendropy
import linkedMSC
import argparse
import random

# replace elements
parser = argparse.ArgumentParser()
parser.add_argument('-t','--tree', default = 'species_tree.tre')
parser.add_argument('-r','--rate', default = 1.26)
parser.add_argument('-N', default = 200)
parser.add_argument('-p', default = 0)
parser.add_argument('-o','--outp', default = 'gene_trees.tre')

args = parser.parse_args()
tree = args.tree
rate = args.rate
N = args.N
p = args.p
outp = args.outp

N= int(N)
p = float(p)
rate = float(rate)

# test 1: 37-taxon species tree
species_tree = dendropy.Tree.get(path = str(tree),schema = "newick")

# --------------------------------------------------------------------------------------------
# generate gene to species map
num_contained = 1
gene_to_species_map = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
        containing_taxon_namespace=species_tree.taxon_namespace,
        num_contained=num_contained)

# --------------------------------------------------------------------------------------------
#  write out the gene to species map for ASTRAL use
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
# gene tree generation
recomb = rate
x = round(N*p)
gene_trees = dendropy.TreeList()

# generate the 0th independent gene tree
gene_tree = dendropy.model.coalescent.contained_coalescent_tree(containing_tree=species_tree,
                                                                 gene_to_containing_taxon_map=gene_to_species_map)
gene_tree.is_rooted = False
gene_trees.append(gene_tree)

# sample x numbers from 1 to N-1 without replacement and sort out from smallest to largest
# which gives the indices of gene trees generated under the MSC model except the 0th gene tree

if x > 0:
    numbers = random.sample(range(1,N), x)
    numbers.sort()
    numbers.append(float('inf'))
    next_indpt = numbers.pop(0)
else:
    next_indpt = -1

for i in range(1,N):
    if i == next_indpt:
        gene_tree = dendropy.model.coalescent.contained_coalescent_tree(containing_tree=species_tree,
                                                                 gene_to_containing_taxon_map=gene_to_species_map)
        gene_tree.is_rooted = False
        gene_trees.append(gene_tree)
        next_indpt = numbers.pop(0)
    else:
        guide_tree = gene_tree
        gene_tree = linkedMSC.multispecies_linked_coalescent(species_tree, guide_tree, recombination_rate=recomb, pop_size = 1)
        gene_tree.is_rooted = False
        gene_trees.append(gene_tree)

gene_trees.write(path=outp,schema='newick')
