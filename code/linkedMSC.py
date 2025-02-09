import dendropy
from dendropy import Tree, TaxonNamespace, Node, Taxon
import numpy as np
from scipy import stats
from scipy.special import comb
import random

def multispecies_linked_coalescent(species_tree, guide_tree, recombination_rate, pop_size = 1, rng = None):
    gene_tree = Tree()
    pop_node_genes_guide = guide_tree.pop_node_genes # a dictionary telling the individuals in each population
    
    height_guiTree = guide_tree.seed_node.distance_from_tip()
    height_spTree = species_tree.seed_node.distance_from_tip()
    gap = height_guiTree - height_spTree # How much taller the guide tree than the species tree is
    
    # pop_node_genes_three is a dictionary that tells what gene lineages reside in each species branch
    # for guide tree lineages and uncoalesced lineages, it is a list
    # for coalesced lineages, it is a dictionary indicates what guide tree lineages they coalesced with
    pop_node_genes_three = {}
    for nd in pop_node_genes_guide:
        pop_node_genes_three[nd] = {'guide_nodes': pop_node_genes_guide[nd],
                                  'coal_nodes': {},
                                  'uncoal_nodes': []}
        
    # lineages in the new locus starts with coalescing with the guide tree
    for nd in pop_node_genes_three:
        if nd.taxon: # a terminal speices branch
            for node in pop_node_genes_three[nd]['guide_nodes']:
                new_node = Node(edge_length = 0.0, taxon = node.taxon)
                pop_node_genes_three[nd]['coal_nodes'][new_node] = node
                
                
    # initialize pop_node_genes for the new gene tree: key is the species branch, value is a list of gene lineages
    # internal branches has no gene lineages for now
    pop_node_genes = {}
    for nd in pop_node_genes_guide:
        coal_nodes = pop_node_genes_three[nd]['coal_nodes']
        pop_node_genes[nd] = []
        for coal_node in coal_nodes:
             pop_node_genes[nd].append(coal_node)
                
    # Start evolving the new gene tree
    for nd in species_tree.postorder_node_iter():
        # root node
        if nd.parent_node is None:
            # the number of lineages in the new locus is more than 1
            if len(pop_node_genes_three[nd]['coal_nodes']) + len(pop_node_genes_three[nd]['uncoal_nodes']) > 1:
                coal_nodes, uncoal_nodes, guide_nodes = linked_coalescent_nodes(pop_node_genes_three[nd]['guide_nodes'],
                                                                   pop_node_genes_three[nd]['coal_nodes'],
                                                                   pop_node_genes_three[nd]['uncoal_nodes'],
                                                                   recombination_rate = recombination_rate,
                                                                   pop_size = pop_size,
                                                                   period = None,
                                                                   now = nd.distance_from_root() + gap,
                                                                   rng = None)
            # only one lineage left in the new locus
            else:
                coal_nodes = pop_node_genes_three[nd]['coal_nodes']
                uncoal_nodes = pop_node_genes_three[nd]['uncoal_nodes']
                
            # set the seed node which is either or not coalesced with the guide tree
            if uncoal_nodes != []:
                gene_tree.seed_node = uncoal_nodes[0]
            else:
                gene_tree.seed_node = list(coal_nodes.keys())[0]
                
        # unroot node
        else:
            coal_nodes, uncoal_nodes, guide_nodes = linked_coalescent_nodes(pop_node_genes_three[nd]['guide_nodes'],
                                                               pop_node_genes_three[nd]['coal_nodes'],
                                                               pop_node_genes_three[nd]['uncoal_nodes'],
                                                               recombination_rate = recombination_rate,
                                                               pop_size = pop_size,
                                                               period = nd.edge_length, # species branch length
                                                               now = nd.distance_from_root() + gap,
                                                               rng = None)
                                
            # update pop_node_genes_three for the parent species branch
            pop_node_genes_three[nd.parent_node]['coal_nodes'].update(coal_nodes)
            pop_node_genes_three[nd.parent_node]['uncoal_nodes'].extend(uncoal_nodes)

            # update pop_node_genes for the parent species branch
            # the resulting gene lineages (coalesecd or uncoalesced with the guide tree) are the input of the parent species branch
            for node in coal_nodes:
                    pop_node_genes[nd.parent_node].append(node)
            for node in uncoal_nodes:
                    pop_node_genes[nd.parent_node].append(node)
                    
    gene_tree.pop_node_genes = pop_node_genes

    return gene_tree



def linked_coalescent_nodes(guide_tree_nodes, coalesced_nodes, uncoalesced_nodes, recombination_rate, pop_size=1, period=None, now=0, rng = None):
    #kg is number of guide tree lineages (coalesced or uncoalesced with a lineage in the new locus)
    #kc is number of lineages coalesced with guide tree (kc <= kg)
    #ku is number of lineages uncoalesced with guide tree
    kg = len(guide_tree_nodes)
    kc = len(coalesced_nodes)
    ku = len(uncoalesced_nodes)
    
    time_remaining = period
        
    while True:
        # stop if it is in the root species branch and there is only one lineage in the new locus
        if time_remaining is None and kc + ku == 1:
            break
        
        #recombination
        if kc >= 1:
            if recombination_rate == 0:
                time_cu = float('inf') # recombination is impossible
            else:
                time_cu = np.random.exponential(scale=1/(recombination_rate * kc)) * pop_size
        else:
            time_cu = float('inf')
        # The waiting time to the next recombination event (exponential distribution)

        #two uncoalesced lineages coalesce
        if ku >= 2:
            time_2u = np.random.exponential(scale=1/(comb(ku, 2))) * pop_size
        else:
            time_2u = float('inf')

        #one uncoalesced lineage coalesces with the guide tree (including coalesced lineage)
        if ku >= 1:
            time_uc = np.random.exponential(scale=1/(ku * kg)) * pop_size
        else:
            time_uc = float('inf')
            
        #two guide tree lineages coalesce (time fixed by tree) (including coalesced lineages)
        if kg >= 2:
            times_coal = {}
            for node in guide_tree_nodes:
                times_coal[node.parent_node] = now - node.parent_node.distance_from_root()
            time_2g = min(times_coal.values())
            coalesced_parent = min(times_coal, key = lambda x:times_coal[x])
            coalesced_children = coalesced_parent.child_nodes()
        else:
            time_2g = float('inf')

        #two coalesced lineages coalesce (time fixed by tree) - subset of '2g' event
        if kc >= 2:
            num_coalescedNode = 0
            for node in coalesced_nodes:
                if coalesced_nodes[node] == coalesced_children[0] or coalesced_nodes[node] == coalesced_children[1]:
                    num_coalescedNode = num_coalescedNode + 1
            if num_coalescedNode == 2:
                time_2c = time_2g
                # only if two candidate nodes for coalescing are in the set of coalesced_nodes, coalescing event happens.
            else:
                time_2c = float('inf')
        else:
            time_2c = float('inf')
        
        events_time = {'cu': time_cu, '2u': time_2u, 'uc': time_uc, '2g': time_2g}
        shortest_time = min(events_time.values())
        
        now = now - shortest_time
        
        event = min(events_time, key = lambda x:events_time[x])
        if time_2c != float('inf') and event == '2g':
            event = '2c'

        if time_remaining is None or shortest_time <= time_remaining:
            #update edge lengths - nodes in new locus only
            for node in uncoalesced_nodes:
                node.edge_length += shortest_time
            for node in coalesced_nodes:
                node.edge_length += shortest_time

            #two coalesced lineages coalesce - coalescence in new locus
            if event == '2c':
                coalesced_pair = []
                
                for node in coalesced_nodes:
                    if coalesced_nodes[node] == coalesced_children[0]:
                        coalesced_pair.append(node)
                    elif coalesced_nodes[node] == coalesced_children[1]:
                        coalesced_pair.append(node)
                # find out which two lineages to coalesce in the set of coalesced_nodes

                coalesced_new = Node(edge_length = 0.0)
                coalesced_new.add_child(coalesced_pair[0])
                coalesced_new.add_child(coalesced_pair[1])
                # Add two lineages to their parent

                coalesced_nodes[coalesced_new] = coalesced_parent
                del coalesced_nodes[coalesced_pair[0]]
                del coalesced_nodes[coalesced_pair[1]]
                # Add coalesced parent node into dic and remove coalesced child nodes in the coalesced lineages
                
                guide_tree_nodes.append(coalesced_parent)
                guide_tree_nodes.remove(coalesced_children[0])
                guide_tree_nodes.remove(coalesced_children[1])
                # coalesce nodes in the guide tree

                kc -= 1
                kg -= 1
                # coalesced_nodes is a dictionary where keys is gene tree nodes and values are guide tree nodes

            #elif shortest_time == time_2g:
            #two guide tree lineages coalesce - no coalescence in new locus
            elif event == '2g':
                for node in coalesced_nodes:
                    if coalesced_nodes[node] == coalesced_children[0] or coalesced_nodes[node] == coalesced_children[1]:
                        coalesced_nodes[node] = coalesced_parent
                # update the mapping of coalesced nodes to the guide tree nodes
                               
                guide_tree_nodes.append(coalesced_parent)
                guide_tree_nodes.remove(coalesced_children[0])
                guide_tree_nodes.remove(coalesced_children[1])
                # coalesce nodes in the guide tree
                kg -= 1

            #elif shortest_time == time_cu:
            #recombination event
            elif event == 'cu':
                candidates = []
                for each in coalesced_nodes:
                    candidates.append(each)
                coalesced_node = random.sample(candidates, 1)[0]

                uncoalesced_nodes.append(coalesced_node)
                del coalesced_nodes[coalesced_node]

                ku += 1
                kc -= 1

            #elif shortest_time == time_2u:
            #two uncoalesced lineages coalesce - coalescence in new locus
            elif event == '2u':
                pair = random.sample(uncoalesced_nodes, 2)
                parent = Node(edge_length = 0.0)
                parent.add_child(pair[0])
                parent.add_child(pair[1])
                uncoalesced_nodes.remove(pair[0])
                uncoalesced_nodes.remove(pair[1])
                uncoalesced_nodes.append(parent)

                ku -= 1
                
            else: # shortest_time = time_uc
                #one uncoalesced lineage coalesces with one guide tree lineage
                isGuideTreeNodeWithCoalescedNodes = False

                uncoalesced_node = random.sample(uncoalesced_nodes, 1)[0]
                guide_tree_node = random.sample(guide_tree_nodes, 1)[0]

                for node in coalesced_nodes:
                    if coalesced_nodes[node] == guide_tree_node:
                        #uncoalesced lineage coalesces with coalesced lineage - coalescence in new locus
                        isGuideTreeNodeWithCoalescedNodes = True

                        child0 = node
                        child1 = uncoalesced_node

                        new_node = Node(edge_length = 0.0)
                        new_node.add_child(child0)
                        new_node.add_child(child1)

                        coalesced_nodes[new_node] = guide_tree_node
                        del coalesced_nodes[child0]
                        uncoalesced_nodes.remove(uncoalesced_node)

                        ku -= 1

                        break
                # uncoalesced lineage coalesce with guide tree lineage which already has lineage coalesced with it.

                #uncoalesced lineage coalesces with guide tree only - no coalescence in new locus
                if not isGuideTreeNodeWithCoalescedNodes:
                    coalesced_nodes[uncoalesced_node] = guide_tree_node
                    uncoalesced_nodes.remove(uncoalesced_node)

                    ku -= 1
                    kc += 1
            
            # adjust the time_remaining left to coalesce
            if time_remaining is not None:
                time_remaining -= shortest_time
            
        else:
            # The next event takes place after the period constraint
            break
        
    if time_remaining is not None and time_remaining > 0:
        for node in uncoalesced_nodes:
            node.edge.length += time_remaining
        for node in coalesced_nodes:
            node.edge.length += time_remaining

    return coalesced_nodes, uncoalesced_nodes, guide_tree_nodes
