import dendropy
from dendropy import Tree, TaxonNamespace, Node, Taxon
import numpy as np
from scipy import stats
from scipy.special import comb
import random

def multispecies_linked_coalescent(species_tree, guide_tree, recombination_rate, pop_size = 1):
    gene_tree = Tree()
    pop_node_genes_guide = guide_tree.pop_node_genes
    
    height_guiTree = guide_tree.seed_node.distance_from_tip()
    height_spTree = species_tree.seed_node.distance_from_tip()
    gap = height_guiTree - height_spTree
    
    # Initialize pop_node_genes_three for the processing
    pop_node_genes_three = {}
    for nd in pop_node_genes_guide:
        pop_node_genes_three[nd] = {'guide_nodes': pop_node_genes_guide[nd],
                                  'coal_nodes': {},
                                  'uncoal_nodes': []}
    for nd in pop_node_genes_three:
        if nd.taxon:
            for node in pop_node_genes_three[nd]['guide_nodes']:
                new_node = Node(edge_length = 0.0, taxon = node.taxon)
                pop_node_genes_three[nd]['coal_nodes'][new_node] = node
                
    # Initialize pop_node_genes for the new gene tree
    pop_node_genes = {}
    for nd in pop_node_genes_guide:
        coal_nodes = pop_node_genes_three[nd]['coal_nodes']
        pop_node_genes[nd] = []
        for coal_node in coal_nodes:
             pop_node_genes[nd].append(coal_node)
                
    # Start growing the new gene tree
    for nd in species_tree.postorder_node_iter():
        # For root node
        if nd.parent_node is None:
            if len(pop_node_genes_three[nd]['coal_nodes']) + len(pop_node_genes_three[nd]['uncoal_nodes']) > 1:
                coal_nodes, uncoal_nodes, guide_nodes = linked_coalescent_nodes(pop_node_genes_three[nd]['guide_nodes'],
                                                                   pop_node_genes_three[nd]['coal_nodes'],
                                                                   pop_node_genes_three[nd]['uncoal_nodes'],
                                                                   recombination_rate = recombination_rate, 
                                                                   pop_size = pop_size, 
                                                                   period = None,
                                                                    now = nd.distance_from_root() + gap)          
            else:
                coal_nodes = pop_node_genes_three[nd]['coal_nodes']
                uncoal_nodes = pop_node_genes_three[nd]['uncoal_nodes']

            if uncoal_nodes != []:
                gene_tree.seed_node = uncoal_nodes[0]
            else:
                gene_tree.seed_node = list(coal_nodes.keys())[0]
                
        # For unroot node
        else:
            coal_nodes, uncoal_nodes, guide_nodes = linked_coalescent_nodes(pop_node_genes_three[nd]['guide_nodes'],
                                                               pop_node_genes_three[nd]['coal_nodes'],
                                                               pop_node_genes_three[nd]['uncoal_nodes'],
                                                               recombination_rate = recombination_rate, 
                                                               pop_size = pop_size, 
                                                               period = nd.edge_length,
                                                                now = nd.distance_from_root() + gap)
            # update pop_node_genes_three for later processing            
            pop_node_genes_three[nd.parent_node]['coal_nodes'].update(coal_nodes)
            pop_node_genes_three[nd.parent_node]['uncoal_nodes'].extend(uncoal_nodes)

            # update pop_node_genes
            for node in coal_nodes:
                    pop_node_genes[nd.parent_node].append(node)
            for node in uncoal_nodes:
                    pop_node_genes[nd.parent_node].append(node)
                    
    gene_tree.pop_node_genes = pop_node_genes
    return gene_tree



def linked_coalescent_nodes(guide_tree_nodes, coalesced_nodes, uncoalesced_nodes, recombination_rate, pop_size=1, period=None, now=0):
    kg = len(guide_tree_nodes)
    kc = len(coalesced_nodes)
    ku = len(uncoalesced_nodes)
    
    time_remaining = period
        
    while kc + ku > 1:

        if kc >= 1:
            if recombination_rate == 0:
                time_cu = float('inf')
            else:
                time_cu = np.random.exponential(scale=1/(recombination_rate * kc)) * pop_size
        else:
            time_cu = float('inf')
        # The waiting time to the next recombination event (exponential distribution)

        if ku >= 2:
            time_2u = np.random.exponential(scale=1/(comb(ku, 2))) * pop_size
        else:
            time_2u = float('inf')

        if ku >= 1:
            time_uc = np.random.exponential(scale=1/(ku * kg)) * pop_size
        else:
            time_uc = float('inf')
            
        if kg >= 2:
            times_coal = {}
            for node in guide_tree_nodes:
                times_coal[node.parent_node] = now - node.parent_node.distance_from_root()
            time_2g = min(times_coal.values())
            coalesced_parent = min(times_coal, key = lambda x:times_coal[x])
            coalesced_children = coalesced_parent.child_nodes()          
        else:
            time_2g = float('inf')

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
            for node in uncoalesced_nodes:
                node.edge_length = node.edge_length + shortest_time
            for node in coalesced_nodes:
                node.edge_length = node.edge_length + shortest_time

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

                kc = kc - 1
                kg = kg - 1
                # coalesced_nodes is a dictionary where keys is gene tree nodes and values are guide tree nodes

            #elif shortest_time == time_2g:
            elif event == '2g':
                for node in coalesced_nodes:
                    if coalesced_nodes[node] == coalesced_children[0] or coalesced_nodes[node] == coalesced_children[1]:
                        coalesced_nodes[node] = coalesced_parent
                # update the mapping of coalesced nodes to the guide tree nodes
                               
                guide_tree_nodes.append(coalesced_parent)
                guide_tree_nodes.remove(coalesced_children[0])
                guide_tree_nodes.remove(coalesced_children[1])
                # coalesce nodes in the guide tree
                kg = kg - 1   

            #elif shortest_time == time_cu:
            elif event == 'cu':
                candidates = []
                for each in coalesced_nodes:
                    candidates.append(each)
                coalesced_node = random.sample(candidates, 1)[0]

                uncoalesced_nodes.append(coalesced_node)
                del coalesced_nodes[coalesced_node]

                ku = ku + 1
                kc = kc - 1

            #elif shortest_time == time_2u:
            elif event == '2u':
                pair = random.sample(uncoalesced_nodes, 2)
                parent = Node(edge_length = 0.0)
                parent.add_child(pair[0])
                parent.add_child(pair[1])
                uncoalesced_nodes.remove(pair[0])
                uncoalesced_nodes.remove(pair[1])
                uncoalesced_nodes.append(parent)

                ku = ku - 1
                
            else: # shortest_time = time_uc
                isGuideTreeNodeWithCoalescedNodes = False

                uncoalesced_node = random.sample(uncoalesced_nodes, 1)[0]
                guide_tree_node = random.sample(guide_tree_nodes, 1)[0]

                for node in coalesced_nodes:
                    if coalesced_nodes[node] == guide_tree_node:
                        isGuideTreeNodeWithCoalescedNodes = True

                        child0 = node
                        child1 = uncoalesced_node

                        new_node = Node(edge_length = 0.0)
                        new_node.add_child(child0)
                        new_node.add_child(child1)

                        coalesced_nodes[new_node] = guide_tree_node
                        del coalesced_nodes[child0]
                        uncoalesced_nodes.remove(uncoalesced_node)

                        ku = ku - 1
                        break
                # uncoalesced lineage coalesce with guide tree lineage which already has lineage coalesced with it.

                if not isGuideTreeNodeWithCoalescedNodes:         
                    coalesced_nodes[uncoalesced_node] = guide_tree_node
                    uncoalesced_nodes.remove(uncoalesced_node)

                    ku = ku - 1
                    kc = kc + 1
            
            # adjust the time_remaining left to coalesce
            if time_remaining is not None:
                time_remaining = time_remaining - shortest_time
            
        else:
            # The next event takes place after the period constraint
            break
        
        if kc + ku == 1 and time_remaining is not None:

            while kg > 1:
                times_coal = {}
                for node in guide_tree_nodes:
                    times_coal[node.parent_node] = now - node.parent_node.distance_from_root()
                time_2g = min(times_coal.values())
                
                if time_2g <= time_remaining:
                    for node in uncoalesced_nodes:
                        node.edge_length = node.edge_length + time_2g
                    for node in coalesced_nodes:
                        node.edge_length = node.edge_length + time_2g

                    coalesced_parent = min(times_coal, key = lambda x:times_coal[x])
                    coalesced_children = coalesced_parent.child_nodes()

                    now = now - time_2g
                    time_remaining = time_remaining - time_2g

                    for node in coalesced_nodes:
                        if coalesced_nodes[node] == coalesced_children[0] or coalesced_nodes[node] == coalesced_children[1]:
                            coalesced_nodes[node] = coalesced_parent
                    # update the mapping of coalesced nodes to the guide tree nodes

                    guide_tree_nodes.append(coalesced_parent)
                    guide_tree_nodes.remove(coalesced_children[0])
                    guide_tree_nodes.remove(coalesced_children[1])
                    # coalesce nodes in the guide tree
                    kg = kg - 1
                    
                else:
                    break
                
            
    if time_remaining is not None and time_remaining > 0:
        for node in uncoalesced_nodes:
            node.edge.length = node.edge_length + time_remaining
        for node in coalesced_nodes:
            node.edge.length = node.edge_length + time_remaining

    return coalesced_nodes, uncoalesced_nodes, guide_tree_nodes