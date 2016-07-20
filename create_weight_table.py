import argparse

import pandas as pd
from Bio import Phylo
from tqdm import tqdm

parser = argparse.ArgumentParser(description="Mappings of classIDs and GIs to counts of occurences of that entity in "
                                             "the classification tree's leaves. Requires the .tre files to be in the "
                                             "same directory as this script and the GI2class files in a certain "
                                             "location (see code).")
parser.add_argument("path_to_tree", help="Path to .tre files.")


args = parser.parse_args()
path_to_tree = args.path_to_tree


def create_mapping(gi_to_class_map, classtree, mapping_filepath, classification_name):
    """
    Create a mapping of GI to a score, scoring the gene by its number of mappings to pathways/classifications.
    If a GI maps to many classification entities it will be scored less.
    :param gi_to_class_map: pandas dataframe containing a gi and classification column.
    :param classtree: tree of the classification system.
    :param mapping_filepath: Path to the mapping file.
    :param classification_name: Name of the classification column in the dataframe (not GI)
    """
    class_ids = pd.Series(
        [node.name for node in tqdm(classtree.root.get_terminals(), desc='Counting multiple class occurences')]
                          ).value_counts().to_dict()
    gi_grouped = gi_to_class_map.groupby(classification_name).groups.items()
    count_gi_to_internal = 0
    with open(mapping_filepath + "gi2count_" + classification_name + ".txt", 'w') as gifile, \
            open(mapping_filepath + "classid2count_" + classification_name + ".txt", 'w') as classificationfile:
        for group, classifications in tqdm(gi_grouped, desc='Computing mapping'):
            if str(group) in class_ids:
                classificationfile.write(str(group) + "\t" + str(class_ids[str(group)]) + "\n")
                for gi in classifications:
                    gifile.write(str(gi) + "\t" + str(class_ids[str(group)]) + "\n")
            else:
                count_gi_to_internal += 1
    print("For", classification_name, count_gi_to_internal, 'GIs were mapped to internal class nodes.')


print('Reading CARD mapping file...')
classification_name = 'card'
card = pd.read_csv(path_to_tree + "card/gi2card-June2016X.map", sep='\t', header=None)
card.columns = ['gi', classification_name]
cardtre = Phylo.read("card1.tre", 'newick')
create_mapping(card, cardtre, path_to_tree, classification_name)
card = None

print('Reading KEGG mapping file...')
classification_name = 'kegg'
kegg = pd.read_csv(path_to_tree + "kegg/gi2kegg-June2016X.map", sep='\t', header=None)
kegg.columns = ['gi', classification_name]
keggtre = Phylo.read("kegg1.tre", 'newick')
create_mapping(kegg, keggtre, path_to_tree, classification_name)
kegg = None

print('Reading SEED mapping file...')
classification_name = 'seed'
seed = pd.read_csv(path_to_tree + "seed/gi2seed-May2015X.map", sep='\t', header=None)
seed.columns = ['gi', classification_name]
seedtre = Phylo.read("seed1.tre", 'newick')
create_mapping(seed, seedtre, path_to_tree, classification_name)
seed = None

print('Reading eggNOG mapping file...')
classification_name = 'eggnog'
eggnog = pd.read_csv(path_to_tree + "eggNOG/gi2eggnog-June2016X.map", sep='\t', header=None)
eggnog.columns = ['gi', classification_name]
eggnogtre = Phylo.read("eggnog1.tre", 'newick')
create_mapping(eggnog, eggnogtre, path_to_tree, classification_name)
eggnog = None

print('Reading InterPro mapping file...')
classification_name = 'interpro'
interpro = pd.read_csv(path_to_tree + "interpro2go/gi2interpro-June2016X.map", sep='\t', header=None)
interpro.columns = ['gi', classification_name]
interprotre = Phylo.read("interpro2go1.tre", 'newick')
create_mapping(interpro, interprotre, path_to_tree, classification_name)
interpro = None

