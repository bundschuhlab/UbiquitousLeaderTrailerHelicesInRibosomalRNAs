#!/usr/bin/env python
import os
import _pickle as pickle
from collections import namedtuple


Node = namedtuple('Node', ['taxid', 'rank', 'parent', 'children'])
name_dict = {
    "d": "domain",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus", 
    "s": "species"
}

def extract_rank_name(feature):
    """ scrub tax names as needed
    """
    feature_rank, _ = feature.split("__")

    feature_rank = name_dict[feature_rank]
    feature_name = feature.lower().replace("\n","")

    return feature_name, feature_rank


def get_lineage(line):
    """ GTDB has repeat taxonomy names for some ranks. Need to dedup
    For duplicate names at diff ranks in the same path, assign name
    to the highest taxonomy rank
    """
    lineage = line.split("\t")[1].split(";")

    # need to de-couple repeat taxonomy names
    lineage_out = []
    lineage_taxonomy_names_only = []
    for l in lineage:
        _, name = l.split("__")
        name = name.lower()
        if name not in lineage_taxonomy_names_only:
            lineage_out.append(l)
            lineage_taxonomy_names_only.append(name)

    return lineage_out


class TaxonomyTree():

    def __init__(self, taxdump_dir, force_reload=False):
        """build gtdb taxonomy tree from bact and arch taxtree

        requires bac120_taxonomy.tsv and ar122_taxonomy.tsv' from
        https://data.gtdb.ecogenomic.org/releases/latest/
        
        >>> tree = gtdb_tree.TaxonomyTree('taxdmp')    
        """
        pickle_fn, bact_nodes_fn, arc_nodes_fn = map(lambda fn: os.path.join(taxdump_dir, fn),
                ('tree.pkl', 'bac120_taxonomy.tsv', 'ar53_taxonomy.tsv'))

        if os.path.exists(pickle_fn) and not force_reload:
            with open(pickle_fn, 'rb') as f:
                self.nodes = pickle.load(f)

        else:
            # add bacteria nodes
            with open(bact_nodes_fn) as f:
                self.nodes = {}
                unseen = {}
                for l in f:
                    lineage = get_lineage(l)

                    for i in range(len(lineage) - 1):
                        parent, child = lineage[i], lineage[i + 1]
                        parent_name, parent_rank = extract_rank_name(parent)
                        child_name, child_rank = extract_rank_name(child)

                        if child_name not in self.nodes:
                            self.nodes[child_name] = Node(
                                    taxid=child_name,
                                    rank=child_rank,
                                    parent=parent_name,
                                    children=set())

                        if parent_name not in self.nodes:
                            # some parent nodes are referenced before they come up as children
                            if parent_name not in unseen:
                                unseen[parent_name] = [set((child_name,)), parent_rank]
                            else:
                                unseen[parent_name][0].add(child_name)
                        else:
                            self.nodes[parent_name].children.add(child_name)

            # add archaea nodes
            with open(arc_nodes_fn) as f:
                for l in f:
                    lineage = get_lineage(l)

                    for i in range(len(lineage) - 1):
                        parent, child = lineage[i], lineage[i + 1]

                        parent_name, parent_rank = extract_rank_name(parent)
                        child_name, child_rank = extract_rank_name(child)

                        if child_name not in self.nodes:
                            self.nodes[child_name] = Node(
                                    taxid=child_name,
                                    rank=child_rank,
                                    parent=parent_name,
                                    children=set())

                        if parent_name not in self.nodes:
                            # some parent nodes are referenced before they come up as children
                            if parent_name not in unseen:
                                unseen[parent_name] = [set((child_name,)), parent_rank]
                            else:
                                unseen[parent_name][0].add(child_name)
                        else:
                            self.nodes[parent_name].children.add(child_name)                

            # second pass creating top nodes
            for parent_name in unseen:
                self.nodes[parent_name] = Node(
                        taxid=parent_name,
                        rank=unseen[parent_name][1],
                        parent=parent_name,
                        children=unseen[parent_name][0])

            with open(pickle_fn, 'wb') as f:
                pickle.dump((self.nodes), f)


    def __getitem__(self, taxid):
        """implements `tree[taxid]` aka dict access
        
        """
        if taxid in self.nodes: return self.nodes[taxid]
        else: raise KeyError('Could not find taxid: {}'.format(taxid))


    def __contains__(self, taxid):
        """implements `taxid in tree`"""
        return taxid in self.nodes
    
    def ascend(self, taxid):
        """walk up the tree and return list of nodes encountered"""
        return [self[taxid]] + (
                self.ascend(self[taxid].parent) if taxid in self and taxid is not self[taxid].parent 
                else [])

    def ascend_taxids(self, taxid):
        """walk up the tree and return list of nodes encountered"""
        return [self[taxid].taxid] + (
                self.ascend_taxids(self[taxid].parent) if taxid in self and taxid is not self[taxid].parent 
                else [])

    def descend(self, taxid):
        """walk down the tree and return set of all child nodes"""
        all_children = {taxid}
        node = self[taxid]

        for child in node.children:
            all_children.update(self.descend(child))
        return all_children
