import pyvolve, dendropy, os
import argparse
import random
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument('-a', type = float, default = 1.062409952497)
parser.add_argument('-b', type = float, default = 0.133307705766)
parser.add_argument('-c', type = float, default = 0.195517800882)
parser.add_argument('-d', type = float, default = 0.223514845018)
parser.add_argument('-e', type = float, default = 0.294405416545)

parser.add_argument('-theta', type = float, default = 0.469075709819)
parser.add_argument('-theta1', type = float, default = 0.558949940165)
parser.add_argument('-theta2', type = float, default = 0.488093447144)

# See https://pbil.univ-lyon1.fr/bpp-doc/bpp-phyl/html/classbpp_1_1GTR.html#details
# for explanation of a,b,c,d,e,theta,theta1,theta2
# See https://www.ideals.illinois.edu/items/55771
# for default values

parser.add_argument('-alpha', help="The shape parameter.", type = float, default =0.370209777709)
parser.add_argument('-n', help="The number of categories to use (n>1).", type = int, default =4)

parser.add_argument('--intp', help="Simulated gene trees.", default ='gene_trees.tre')
parser.add_argument('--seqlen', help="Simulated sequence length", default='mix')
parser.add_argument('--outp', help="directory to store simulated sequences", default='./')

parser.add_argument('--seed', type=int, help="Random seed for reproducibility.", default=None)

args = parser.parse_args()

intp = args.intp
seqlen = args.seqlen
outp = args.outp

a = args.a
b = args.b
c = args.c
d = args.d
e = args.e

theta = args.theta
theta1 = args.theta1
theta2 = args.theta2

alpha = args.alpha
n = args.n
seed = args.seed

# Set the seed if provided
if seed is not None:
    random.seed(seed)
    print(f"Random seed set to: {seed}")
else:
    print("No seed provided. Using default random behavior.")

## read the biological gene trees
#biotrees = dendropy.TreeList.get(path=bio, schema="newick")
# Read the simulated gene trees
trees = dendropy.TreeList.get(path=intp, schema="newick")
num = len(trees)


### Rescale branch lengths of simulated gene trees to match biological trees
## the number of branches
#num_taxa = biotrees[0].__len__()
#num_edges = 2*num_taxa - 3
#
## calculate the mean of the average branch lengths of 424 biological gene trees
#aveBranchLengths = []
#for t in biotrees:
#    sumBLs = t.length()
#    aveBranchLengths.append(sumBLs/num_edges)
#aveBLmean = np.mean(aveBranchLengths)
#
## Scale branch lengths for each simulated tree
#for t in trees:
#    aBL = t.length()/num_edges # calculate the average branch length of a true gene tree
#    scale = aveBLmean/aBL
#    for node in t.__iter__():
#        if node.edge.length is not None:
#            node.edge.length = node.edge.length * scale
        
###########################################
## Evolve sequences using scaled gene trees

# Split indices of gene trees randomly into two halfs
def split_indices(num):
    indices = list(range(num))
    random.shuffle(indices)
    
    mid = num // 2
    first_half = indices[:mid]
    
    # If n is even, the first half will have one less index than the second half
    return first_half

# Construct nucleotide model
def model_gen(a,b,c,d,e,theta,theta1,theta2):
    # Frequencies, in order A C G T.
    piA = theta1*(1-theta)
    piC = (1-theta2)*theta
    piG = theta2*theta
    piT = (1-theta1)*(1-theta)
    freqs = [piA, piC, piG, piT]
    
    # Mutation rates
    P = 2*(a*piC*piT
           +b*piA*piT
           +c*piG*piT
           +d*piA*piC
           +e*piC*piG
           +piA*piG)
    AC = d*piC/P
    AG = piG/P
    AT = b*piT/P
    CG = e*piG/P
    CT = a*piT/P
    GT = c*piT/P
    mu = {"AC":AC, "AG":AG, "AT":AT, "CG":CG, "CT":CT, "GT":GT}
    
    # Construct nucleotide model with mutation rates and frequencies.
    nuc_model = pyvolve.Model("nucleotide", {"mu":mu, "state_freqs":freqs}, alpha=alpha, num_categories=n)
    return nuc_model

nuc_model = model_gen(a,b,c,d,e,theta,theta1,theta2)

# Half of gene trees evolve sequence of length 500
for i in range(num):
    if seed is not None:
        random.seed(seed)
#    trees[i].is_rooted = None
#    trees[i].write(path='tmp_tree.tre',schema='newick')
#    phylogeny = pyvolve.read_tree(file = 'tmp_tree.tre')
    # convert to newick string
    newick_str = trees[i].as_string(schema="newick")

    # read tree directly from string
    phylogeny = pyvolve.read_tree(tree=newick_str)
    # In nucleotide models, branch lengths represent mean number of subsititutions per unit time
    
    if seqlen == "mix":
        first_half = split_indices(num)
#        print('Indices of the first half of gene trees:', first_half)
        if i in set(first_half):
            my_partition = pyvolve.Partition(models = nuc_model, size = 500)
        else:
            my_partition = pyvolve.Partition(models = nuc_model, size = 1000)
    else:
        my_partition = pyvolve.Partition(models = nuc_model, size = int(seqlen))
    # size is the number of sites to evolve in this partition

    my_evolver = pyvolve.Evolver(tree = phylogeny, partitions = my_partition)
    # partition =[partition1. partition2, partition3] can evolve several partitions

    # specify the path for the directory – make sure to surround it with quotation marks
    # ./ stands for the current working directory
    path = f'./{outp}/{i}'

    # create new single directory for the convenience of the following gene tree estimation
    os.makedirs(path, exist_ok=True)

    my_evolver(seqfile = f"{path}/seq{i}.fasta", ratefile = None, infofile = None, seed = seed) # , seqfmt = "phylip"
