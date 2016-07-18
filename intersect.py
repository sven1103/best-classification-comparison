import argparse

import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get GI's intersecting with the different classification systems:"
                                                 "EGGnog, KEGG, SEED, InterPro, CARD).")
    parser.add_argument("seedmap", help="Path to SEED map file.")
    parser.add_argument("keggmap", help="Path to KEGG map file.")
    parser.add_argument("eggNOGmap", help="Path to eggNOG map file.")
    parser.add_argument("cardmap", help="Path to CARD map file.")
    parser.add_argument("interpromap", help="Path to Interpro map file.")

    args = parser.parse_args()

    print('Reading KEGG mapping file...')
    kegg = pd.read_csv(args.keggmap, sep='\t')
    kegg.columns = ['gi', 'kegg']
    print('Reading SEED mapping file...')
    seed = pd.read_csv(args.seedmap, sep='\t')
    seed.columns = ['gi', 'seed']
    print('Reading eggNOG mapping file...')
    eggnog = pd.read_csv(args.eggNOGmap, sep='\t')
    eggnog.columns = ['gi', 'eggnog']
    print('Reading CARD mapping file...')
    card = pd.read_csv(args.cardmap, sep='\t')
    card.columns = ['gi', 'card']
    print('Reading InterPro mapping file...')
    interpro = pd.read_csv(args.interpromap, sep='\t')
    interpro.columns = ['gi', 'interpro']

    keggseed = pd.merge(kegg, seed, how='inner', on='gi')
    keggnog = pd.merge(kegg, eggnog, how='inner', on='gi')
    keggcard = pd.merge(kegg, card, how='inner', on='gi')
    kegginterpro = pd.merge(kegg, interpro, how='inner', on='gi')
    seedeggnog = pd.merge(seed, eggnog, how='inner', on='gi')
    seedcard = pd.merge(seed, card, how='inner', on='gi')
    seedinterpro = pd.merge(seed, interpro, how='inner', on='gi')
    eggnogcard = pd.merge(eggnog, card, how='inner', on='gi')
    eggnoginterpro = pd.merge(eggnog, interpro, how='inner', on='gi')
    cardinterpro = pd.merge(card, interpro, how='inner', on='gi')

