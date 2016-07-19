
# coding: utf-8

# In[2]:

get_ipython().magic('load_ext Cython')
get_ipython().magic('matplotlib inline')
from Bio import Phylo
import os
import matplotlib.pylab as plt
import pandas as pd
import collections
from sklearn.metrics import adjusted_rand_score,jaccard_similarity_score
import numpy as np
import warnings

PATH_TO_TREE = '/home/sven1103/megan/class/resources/files/'

PATH_TO_MAP = '/home/sven1103/ownCloud/notebooks/BEST/'

PHYLO_TYPE = "newick"


def get_trees():
    trees = {}
    for file_name in os.listdir(PATH_TO_TREE):
        if file_name.split(".")[1] in "tre":
            print(file_name)
            tree_path = PATH_TO_TREE + file_name
            print(tree_path)
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
        
        # All nodes in all levels
        depth_series.groupby(depth_series).size()
        
        # Node sum
        depth_series.groupby(depth_series).size().sum()
        
        # Print all stats
        print(str.format("{0}: {1}", "Number of leaves", number_terminals))
        print(str.format("{0}\n {1}", "Levels and number of nodes", depth_series.groupby(depth_series).size()))
        print(str.format("{0}: {1}", "Internal nodes", depth_series.groupby(depth_series).size().sum() - number_terminals))
        print("--------------------------------------")
        
        leave_list = list(map(lambda x: int(x.name), dict_trees[key].get_terminals()))
        leave_series = pd.Series(leave_list)
        leave_series.name = "leaves"
        leave_counting_dict = dict(collections.Counter(leave_series.groupby(leave_series).size()))
        
        plt.bar(leave_counting_dict.keys(), leave_counting_dict.values())
        plt.yscale("log")
        plt.title(key.split(".")[0])
        plt.xlabel("Number of pathway assignments")
        plt.ylabel("Number of gene families")
        plt.savefig(PATH_TO_MAP+key.split(".")[0]+".pdf")
        plt.close()


# In[3]:

# Read in the trees
dict_trees = get_trees()


# In[3]:

get_basic_stats(dict_trees)


# In[4]:

leave_list = list(map(lambda x: int(x.name), dict_trees["seed1.tre"].get_terminals()))


# ## Complemetary classifications

# In[4]:

leaves_kegg = dict_trees["kegg1.tre"].get_terminals()
leaves_parents = dict_trees["kegg1.tre"].get_nonterminals()
subnodes_kegg = {}

for parent_node in leaves_parents:
    #print(parent_node.confidence)
    if parent_node.confidence not in subnodes_kegg:
        subnodes_kegg[parent_node.confidence] = []
    for leave in leaves_kegg:
        if leave in parent_node:
            subnodes_kegg[parent_node.confidence].append(leave.name)


# In[5]:

leaves_seed = dict_trees["seed1.tre"].get_terminals()
leaves_parents = dict_trees["seed1.tre"].get_nonterminals()
subnodes_seed = {}

for parent_node in leaves_parents:
    #print(parent_node.confidence)
    if parent_node.confidence not in subnodes_seed:
        subnodes_seed[parent_node.confidence] = []
    for leave in leaves_seed:
        if leave in parent_node:
            subnodes_seed[parent_node.confidence].append(leave.name)


# ### Vorsicht!
# leere Listen entfernen!

# In[ ]:




# In[7]:

PATH_TO_MAP


# In[6]:

kegg_map = pd.read_table(PATH_TO_MAP + "gi2kegg-June2016X.map", sep="\t", names=["gi","id"])


# In[7]:

seed_map = pd.read_table(PATH_TO_MAP + "gi2seed-May2015X.map", sep="\t", names=["gi","id"])


# In[10]:

len(kegg_map.loc[kegg_map["id"]==1,:])


# In[8]:

def get_dict_map(mapping_file):
    """Converts a mapping file with gi number and class identifier
    to a dictionary with: {class_id:[gi_ids]}
    """
    
    dict_map = {}
    for key, df in mapping_file.groupby("id"):
        dict_map[key] = df["gi"].tolist()
    return dict_map


# In[9]:

dict_map_kegg = get_dict_map(kegg_map)


# In[10]:

dict_map_seed = get_dict_map(seed_map)


# In[22]:

get_ipython().run_cell_magic('cython', '', 'cdef list convert_to_gi_list(list integer_list, dict dict_map):\n    tmp = []\n    if not integer_list:\n        return tmp\n    for integer in integer_list:\n        if int(integer) in dict_map:\n            tmp.extend(dict_map[int(integer)])\n    return tmp\n\ncdef double weighted_jaccard_index(list list_1, list list_2):\n    cdef double index\n    if not list_1 or not list_2:\n        return 0\n    set_1 = set(list_1)\n    set_2 = set(list_2)\n    joint_set = set_1.intersection(set_2)\n    index = len(joint_set)/(min(len(set_1), len(set_2))+0.0)\n    return index\n\n\ndef compare_classifications(class_sets_1, class_sets_2):\n    cdef int counter\n    cdef double total_size, total_inner, max_jaccard\n    cdef list jaccard_index_list\n    jaccard_index_list = []\n    counter=1\n    total_size = len(class_sets_1)\n    for leave_list in class_sets_1:\n        max_jaccard = 0\n        for leave_list_other in class_sets_2:\n            current_jaccard = weighted_jaccard_index(leave_list,leave_list_other)\n            max_jaccard = max(max_jaccard,current_jaccard)\n        jaccard_index_list.append(max_jaccard)\n        print(str.format("{0}/{1} compared.", counter,total_size))\n        counter += 1\n    return jaccard_index_list\n\ncdef list make_set_of_gi(subnodes, dict dict_map):\n    set_list = []\n    for leave_list in subnodes.values():\n        set_list.append(convert_to_gi_list(leave_list, dict_map))\n    return set_list')


# In[23]:

kegg_sets = make_set_of_gi(subnodes_kegg, dict_map_kegg)


# In[24]:

seed_sets = make_set_of_gi(subnodes_seed, dict_map_seed)


# In[16]:

jaccard_list = compare_classifications(kegg_sets, seed_sets)


# In[18]:

jaccard_list

