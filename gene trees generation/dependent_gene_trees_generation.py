import dendropy
import linkedMSC
import argparse

# replace elements
parser = argparse.ArgumentParser()
parser.add_argument('-t','--tree', default = 'species_tree.tre')
parser.add_argument('-r','--rate', default = 1e6)
parser.add_argument('-n','--num', default = 200)
parser.add_argument('-o','--outp', default = 'gene_trees.tre')

args = parser.parse_args()
tree = args.tree
rate = args.rate
num = args.num
outp = args.outp

# test 1: 37-taxon species tree
species_tree = dendropy.Tree.get(path = str(tree),schema = "newick")

# --------------------------------------------------------------------------------------------
# generate gene to species map
num_contained = 1
gene_to_species_map = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
        containing_taxon_namespace=species_tree.taxon_namespace,
        num_contained=num_contained)

# --------------------------------------------------------------------------------------------
# write out the gene to species map for ASTRAL use
#f = open("gene_to_species_map.txt","w")
#
#species = list(species_tree.taxon_namespace)
#individuals = list(gene_to_species_map)
#
#num_species = len(species)
#individuals_start_index = 0
#
#if type(num_contained) == int:
#    num_contained = [num_contained] * num_species
#
#for i in range(num_species):
#    f.write(species[i].label + ':')
#    num_individuals = num_contained[i]
#    if num_individuals == 1:
#        f.write(individuals[individuals_start_index].label.replace(" ", "_"))
#    else:       
#        for j in range(num_individuals - 1):
#            f.write(individuals[individuals_start_index + j].label.replace(" ", "_") + ',')
#        f.write(individuals[individuals_start_index + j + 1].label.replace(" ", "_"))
#    individuals_start_index = individuals_start_index + num_individuals
#    f.write('\n') 
#
#f.close()

# --------------------------------------------------------------------------------------------
# generate gene trees and store in a file
from dendropy import TreeList
gene_trees = TreeList()
guide_tree = dendropy.model.coalescent.contained_coalescent_tree(containing_tree=species_tree,
                                                                    gene_to_containing_taxon_map=gene_to_species_map)
gene_trees.append(guide_tree)
for i in range(int(num)-1):
    gene_tree = linkedMSC.multispecies_linked_coalescent(species_tree, guide_tree, recombination_rate=float(rate), pop_size = 1)
    gene_tree.is_rooted = False
    #gene_tree.update_bipartitions(suppress_unifurcations=False)
    gene_trees.append(gene_tree)
    guide_tree = gene_tree
gene_trees.write(path=outp,schema='newick')
