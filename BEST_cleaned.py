from Bio import Phylo
import networkx as nxd
import os
import matplotlib.pylab as plt
import seaborn as sns
import pandas as pd

PATH_TO_TREE = '/home/sven1103/megan/class/resources/files/'

PHYLO_TYPE = "newick"


def get_trees():
    trees = {}
    for file_name in os.listdir(PATH_TO_TREE):
        if file_name.split(".")[1] in "tre":
            tree_path = PATH_TO_TREE + file_name
            phylo_tree = Phylo.read(tree_path, PHYLO_TYPE)
            trees.update({file_name:phylo_tree})
    return trees


def get_basic_stats(dict_trees):
    for key,_ in dict_trees.items():
        print("Stats for", key)
        number_terminals = dict_trees[key].count_terminals()
        depths_dict = dict_trees[key].depths(unit_branch_lengths=True)
        depths = list(depths_dict.values())
        depth_series = pd.Series(depths)
        depth_series.name = "depth_levels"
        
        # Print all stats
        # Print the number of leaves
        print(str.format("{0}: {1}", "Number of leaves", number_terminals))
        
        # Print the depth information
        print(str.format("{0}\n {1}", "Levels and number of nodes",
                         depth_series.groupby(depth_series).size()))
        
        # Print the internal nodes
        print(str.format("{0}: {1}", "Internal nodes",
                         depth_series.groupby(depth_series).size().sum() - number_terminals))
        
        print("--------------------------------------")


