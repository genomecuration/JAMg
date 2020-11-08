#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import os, sys, re
import logging
import argparse
import collections
import numpy
import time

import TNode

logger = logging.getLogger(__name__)


class TGraph:

 
    def __init__(self, gene_id):

        self.node_cache = dict()
        self.gene_id = gene_id

    
    def get_node(self, transcript_id, loc_node_id, node_seq):
        
        """
        Instantiates Node objects, and stores them in a graph.

        *** use this method for instantiating all Node objects ***

        use clear_node_cache() to clear the graph
                
        """
                
        logger.debug("{}\t{}".format(loc_node_id, node_seq))
        
        if len(node_seq) == 0:
            raise RuntimeError("Error, non-zero length node_seq required for parameter")


        if loc_node_id in self.node_cache:
            node_obj = self.node_cache[ loc_node_id ]
            node_obj.add_transcripts(transcript_id)
            if node_obj.seq != node_seq:
                errmsg = "ERROR: have conflicting node sequences for {} node_id: {}\n".format(self.get_gene_id(),
                                                   loc_node_id) + "{}\n vs. \n{}".format(node_obj.seq, node_seq)
                logger.critical(errmsg)
                
                raise RuntimeError(errmsg)
            else:
                return node_obj
            
        else:
            # instantiate a new one
            node_obj = TNode.TNode(self, transcript_id, loc_node_id, node_seq)
            self.node_cache[ loc_node_id ] = node_obj
            return node_obj



    def get_all_nodes(self):
        return list(self.node_cache.values())
    
    def clear_node_cache(self):
        """
        clears the graph
        """
        self.node_cache.clear()
    
    def clear_touch_settings(self):
        """
        clear the touch settings for each of the nodes
        """

        for node in self.get_all_nodes():
            node.clear_touch()
    


    def add_edges(self, from_nodes_list, to_nodes_list):

        for from_node in from_nodes_list:
            for to_node in to_nodes_list:
                from_node.add_next_node(to_node)
                to_node.add_prev_node(from_node)

    def prune_edges(self, from_nodes_list, to_nodes_list):

        for from_node in from_nodes_list:
            for to_node in to_nodes_list:
                from_node.remove_next_node(to_node)
                to_node.remove_prev_node(from_node)
    

    def prune_node(self, node):
        logger.debug("pruning node: {}".format(node))
        self.prune_edges(node.get_prev_nodes(), [node])
        self.prune_edges([node], node.get_next_nodes())
        node.dead = True
        self.node_cache.pop(node.get_loc_id())
    

    def retrieve_node(self, node_id):
        """
        does not instantiate, only retrieves.
        If loc_node_id is not in the graph, returns None
        """

        if node_id in self.node_cache:
            return self.node_cache[node_id]
        else:
            return None
        
        
    def get_gene_id(self):
        return self.gene_id


    def draw_graph(self, filename):

        logger.debug("drawing graph: {}".format(filename))
        ofh = open(filename, 'w')

        ofh.write("digraph G {\n")

        for node_id in self.node_cache:
            node = self.node_cache[node_id]
            node_seq = node.get_seq()
            gene_node_id = node.get_gene_node_id()
            next_nodes = node.get_next_nodes()
            node_seq_len = len(node_seq)

            max_len_show = 50
            max_len_show_half = int(max_len_show/2)

            if node_seq_len > max_len_show:
                node_seq = node_seq[0:max_len_show_half] + "..." + node_seq[(node_seq_len-max_len_show_half):node_seq_len]
            
            ofh.write("{} [label=\"{}:Len{}:T{}:{}\"]\n".format(node.get_id(), gene_node_id, node_seq_len,
                                                                node.get_topological_order(), node_seq))
            
            for next_node in next_nodes:
                ofh.write("{}->{}\n".format(node.get_id(), next_node.get_id()))

        ofh.write("}\n") # close it

        ofh.close()
    
    def __repr__(self):
        txt = ""
        for node in self.get_all_nodes():
            txt += node.toString() + "\n"

        return txt
