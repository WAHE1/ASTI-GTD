import dendropy
import argparse

# replace elements
parser = argparse.ArgumentParser()
parser.add_argument('-t','--sptree',default = 'species_tree.tre')
parser.add_argument('-i','--iftree',default = 'inferred_tree.tre')
parser.add_argument('-o','--outp',default = 'RFdist.txt')

args = parser.parse_args()
outp = args.outp
sptree = args.sptree
iftree = args.iftree

tns = dendropy.TaxonNamespace() 
species_tree = dendropy.Tree.get(
        path = sptree,
        schema = "newick",
        rooting='force-unrooted',
        taxon_namespace=tns)
inferred_tree = dendropy.Tree.get(
        path = iftree,
        schema = "newick",
        taxon_namespace=tns)

RF_max = (len(tns)-3)*2

RF_distance = dendropy.calculate.treecompare.unweighted_robinson_foulds_distance(species_tree, inferred_tree)
print(f'The RF distance between {sptree} and {iftree} is {RF_distance}')
print(f'The Normalised RF distance between {sptree} and {iftree} is {RF_distance/RF_max}')

f = open(outp,'a')
f.write(f'{RF_distance}\n')
f.close()
