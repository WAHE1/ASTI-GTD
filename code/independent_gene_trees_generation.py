import dendropy
import argparse
import random
parser = argparse.ArgumentParser()
parser.add_argument('-t','--tree', default = 'species_tree.tre')
parser.add_argument('-n','--num', type = int, default = 200)
parser.add_argument('-o','--outp', default = 'gene_trees.tre')
parser.add_argument('--seed', type=int, help="Random seed for reproducibility.", default=None)

args = parser.parse_args()
tree = args.tree
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
species_tree = dendropy.Tree.get(path = tree,schema = "newick")
# generate gene to species map
num_contained = 1
gene_to_species_map = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
        containing_taxon_namespace=species_tree.taxon_namespace,
        num_contained=num_contained)

# generate gene trees and store in a file
gene_trees = dendropy.TreeList()
for i in range(num):
    gene_tree = dendropy.model.coalescent.contained_coalescent_tree(containing_tree=species_tree,
                                                                    gene_to_containing_taxon_map=gene_to_species_map,
                                                                    rng=rng)
    gene_trees.append(gene_tree)
    
## To make taxon lable from taxon_1 to taxon
#for i in range(num):
#    for t in gene_trees[i].leaf_nodes():
#        t.taxon.label = t.taxon.label.split( )[0]
        
gene_trees.write(path=outp,schema='newick')
