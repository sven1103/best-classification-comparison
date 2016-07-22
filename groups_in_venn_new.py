from Bio import Phylo
import pandas as pd

path_to_tree = "../../"
PHYLO_TYPE = "newick"

def all_parents(tree):
    parents = {}
    for clade in tree.find_clades():
        for child in clade:
        	if (child in parents):
        		print("ERROR")
        	parents[child] = clade
    return parents

def lookup_names(tree):
	names = {}
	count = 0
	for clade in tree.find_clades():
		if (clade.name == None):
			if (clade.confidence == None):
				count += 1
			elif (clade.confidence in names):
				names[clade.confidence].append(clade)
			else:
				names[clade.confidence] = [clade]
		elif (clade.name in names):
			names[clade.name].append(clade)
		else:
			names[clade.name] = [clade]
	print(count)
	return names

print('Reading KEGG mapping file...')
kegg = pd.read_csv(path_to_tree+"kegg/gi2kegg-June2016X.map", sep='\t', header=None)
kegg.columns = ['gi', 'kegg']

print('Reading eggNOG mapping file...')
eggnog = pd.read_csv(path_to_tree+"eggNOG/gi2eggnog-June2016X.map", sep='\t', header=None)
eggnog.columns = ['gi', 'eggnog']

intersect = pd.merge(kegg,eggnog,"inner")
left = pd.merge(kegg,eggnog,"left")
right = pd.merge(kegg,eggnog,"right")

k_i = set(intersect["kegg"])
e_i = set(intersect["eggnog"])
l = set(left["kegg"])
r = set(right["eggnog"])

kegg_tree = Phylo.read(path_to_tree + "kegg/kegg.tre", PHYLO_TYPE)
egg_tree = Phylo.read(path_to_tree + "eggNOG/eggnog.tre", PHYLO_TYPE)

for key in k_i:
