# ASTI-GTD

## linkedMSC.py
Simulate a dependent gene tree given a guide tree constrained by a species tree (in coalescent units)

To run linkedMSC.py, package `dendropy` is needed (see https://jeetsukumaran.github.io/DendroPy/ for more details about how to install it). Copy linkedMSC.py to the working directory.
### An example:
```
import linkedMSC
gene_tree = linkedMSC.multispecies_linked_coalescent(species_tree, guide_tree, recombination_rate = 1e6, pop_size = 1, rng = 1234)
```

## dependent_gene_trees_generation.py
Simulate a set of dependent gene trees within a species tree

### Arguments:
`-r [a number]`: recombination rate

`-n [an integer]`: number of gene trees

`-t [a newick tree filename]`: species tree

`-o [a newick tree filename]`: output gene trees
### An example:
Simulate 200 gene trees under the recombination rate 1e6. The constraint species tree is `sptree.tre`, and the simulated gene trees is written in `gene_trees.tre`.
```
python dependent_gene_trees_generation.py -t sptree.tre -r 1e6 -n 200 -o gene_trees.tre --seed None
```

## independent_gene_trees_generation.py
Simulate a set of independent gene trees within a species tree under the MSC model

### Arguments:

`-n [an integer]`: number of gene trees

`-t [a newick tree filename]`: species tree

`-o [a newick tree filename]`: output gene trees

`--seed [an integer]`: random seed
### An example:
Simulate 200 independent gene trees under the MSC model. The constraint species tree is `sptree.tre`, and the simulated gene trees is written in `gene_trees.tre`.
```
python independent_gene_trees_generation.py -t sptree.tre -n 200 -o gene_trees.tre --seed None
```

## fixRbar_inAdpt_geneTrees_generation.py/fixR_inAdpt_geneTrees_generation.py
Simulate multiple independent sets of dependent gene trees by fixing the overall recombination rate between all trees/fixing the recombination rate between trees in each set.

### Arguments:
`-r [a number]`: recombination rate

`-N [an integer]`: number of gene trees

`-t [a newick tree filename]`: species tree

`-o [a newick tree filename]`: output gene trees

`-p [a number between 0 and 1]`: proportion of independent gene trees
### An example:
```
python fixR_inAdpt_genetrees_generation.py -t sptree.tre -r 1e6 -N 200 -p 0 -o gene_trees.tre
```

## interval_sample_trees.py
Sample gene trees that are a certain distance apart from a set of trees.

### Arguments:
`-i [a newick tree filename]`: input trees

`-s [an integer]`: the width of the interval (sample every s trees)

`-o [a newick tree filename]`: output gene trees
### An example:
```
python interval_sample_trees.py -i trees.tre -s 10 -o sample_trees.tre
```

## sequence_generation.py
### Arguments:
`--seqlen [an integer]`: sequence length

`--outp [directory]`: directory to save simulated sequences
### An example:
```
python sequences_generation.py -a 1.062409952497 -b 0.133307705766 -c 0.195517800882 -d 0.223514845018 -e 0.294405416545 -theta 0.469075709819 -theta1 0.558949940165 -theta2 0.488093447144 -alpha 0.370209777709 -n 4 --intp gene_trees.tre --seqlen 500 --outp data
```
See information about `a, b, c, d, e, theta, theta1` and `theta2` from https://pbil.univ-lyon1.fr/bpp-doc/bpp-phyl/html/classbpp_1_1GTR.html#details.

## RF_distance.py
Calculate the Robinson Foulds distance between two trees
### Arguments:
`-t [a newick tree filename]`: a tree

`-i [a newick tree filename]`: a tree

`-o [a txt filename]`: file used to write out the RF distance
### An example:
```
python RF_distance.py -t sptree.tre -i iftree.tre -o RFdist.txt
```


