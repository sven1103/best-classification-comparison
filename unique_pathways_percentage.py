import argparse
from Bio import Phylo
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get the number of unique pathways/classes leaves in the different "
                                                 "classification systems: EGGnog, KEGG, SEED, InterPro, CARD).")
    parser.add_argument("seedtre", help="Path to SEED tre file.")
    parser.add_argument("keggtre", help="Path to KEGG tre file.")
    parser.add_argument("eggnogtre", help="Path to eggNOG tre file.")
    parser.add_argument("cardtre", help="Path to CARD tre file.")
    parser.add_argument("interprotre", help="Path to InterPro tre file.")
    parser.add_argument("interpromap", help="Path to InterPro map file.")

    args = parser.parse_args()

    classifications = {
        'SEED': Phylo.read(args.seedtre, 'newick'),
        'KEGG': Phylo.read(args.keggtre, 'newick'),
        'eggNOG': Phylo.read(args.eggnogtre, 'newick'),
        'CARD': Phylo.read(args.cardtre, 'newick'),
        'InterPro': Phylo.read(args.interprotre, 'newick'),
    }

    def getuniquepercentage(classificationtree):
        totalleafcount = classificationtree.count_terminals()
        leaves = classificationtree.get_terminals()
        s = []
        for clade in leaves:
            if clade.name not in s:
                s.append(clade.name)
        return {
            'totalleafcount': totalleafcount,
            'uniqueleafcount': len(s),
            'percentage': float(len(s)) / float(totalleafcount)
        }

    for classification, tree in classifications.items():
        res = getuniquepercentage(tree)
        print("Results for", classification)
        print("\tTotal leaf count:", res['totalleafcount'])
        print("\tNumber of uniquely occuring leaves:", res['uniqueleafcount'])
        print("\tThe percentage of unique leaves:", res['percentage']*100, "%")
        print("#################################################################")
        print()

    # Have a closer look to the InterPro subtrees 'molecular function', 'cellular components', 'biological process' and
    # 'Unclassified'
    mapping = pd.read_csv(args.interpromap, sep='\t', names=['nodeid', 'name', 'description'])
    subtrees = [subtree for subtree in classifications['InterPro'].root.clades if subtree.confidence is not None]
    #classifications['InterPro'].root.clades[]
    print("Results for InterPro")
    for subtree in subtrees:
        print("\tSubtree: ",  mapping[mapping['nodeid'] == subtree.confidence].name.values[0])
        res = getuniquepercentage(subtree)
        print("\tTotal leaf count:", res['totalleafcount'])
        print("\tNumber of uniquely occuring leaves:", res['uniqueleafcount'])
        print("\tThe percentage of unique leaves:", res['percentage'] * 100, "%")
        print("\t########################################################")

    print("#################################################################")

