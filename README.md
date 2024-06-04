# ASTI-GTD

## linkedMSC.py
Given a guide tree, a gene tree is simulated within a species tree.
### An example:
```
import linkedMSC
gene_tee = linkedMSC.multispecies_linked_coalescent(species_tree, guide_tree, recombination_rate)
```

## dependent_gene_trees_generation.py
### Arguments:
`-r [a number]`: recombination rate

`-n [an integer]`: number of gene trees

`-t [a newick tree filename]`: species tree

`-o [a newick tree filename]`: output gene trees

## fixRbar_indptAdpt_geneTrees_generation.py/fixR_indptAdpt_geneTrees_generation.py
### Arguments:
`-r [a number]`: recombination rate

`-N [an integer]`: number of gene trees

`-t [a newick tree filename]`: species tree

`-o [a newick tree filename]`: output gene trees

`-p [a number between 0 and 1]`: proportion of independent gene trees


