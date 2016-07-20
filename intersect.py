import argparse
import pandas as pd
from matplotlib_venn import venn2
from matplotlib_venn import venn3
from matplotlib import pyplot as plt
from tqdm import tqdm

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get GI's intersecting with the different classification systems:"
                                                 "EGGnog, KEGG, SEED, InterPro, CARD).")
    parser.add_argument("path_to_tree", help="Path to tree files.")
    

    args = parser.parse_args()
    path_to_tree = args.path_to_tree

    print('Reading KEGG mapping file...')
    kegg = pd.read_csv(path_to_tree+"kegg/gi2kegg-June2016X.map", sep='\t', header=None)
    kegg.columns = ['gi', 'kegg']
    kegg = kegg.drop_duplicates(subset='gi')
    print('Reading SEED mapping file...')
    seed = pd.read_csv(path_to_tree+"seed/gi2seed-May2015X.map", sep='\t', header=None)
    seed.columns = ['gi', 'seed']
    seed = seed.drop_duplicates(subset='gi')
    print('Reading eggNOG mapping file...')
    eggnog = pd.read_csv(path_to_tree+"eggNOG/gi2eggnog-June2016X.map", sep='\t', header=None)
    eggnog.columns = ['gi', 'eggnog']
    eggnog = eggnog.drop_duplicates(subset='gi')
    print('Reading CARD mapping file...')
    card = pd.read_csv(path_to_tree+"card/gi2card-June2016X.map", sep='\t', header=None)
    card.columns = ['gi', 'card']
    card = card.drop_duplicates(subset='gi')
    print('Reading InterPro mapping file...')
    interpro = pd.read_csv(path_to_tree+"interpro2go/gi2interpro-June2016X.map", sep='\t', header=None)
    interpro.columns = ['gi', 'interpro']
    interpro = interpro.drop_duplicates(subset='gi')

    names = ["KEGG","SEED","eggNOG","CARD","interpro2go"]
    outer_sets = [kegg,seed,eggnog,card,interpro]

    intersects = [[0 for x in range(len(names))] for y in range(len(names))]
    count = 0
    indexes = [[-1 for x in range(len(names))] for y in range(len(names))]
    inters = []
    candice = [[0 for x in range(len(names))]for y in range(len(names))]
    for i in tqdm(range(0,len(names)-1)):
        for j in range(i+1,len(names)):
            print("Running merge on %s with %s..." % (names[i], names[j]))
            inter = pd.merge(outer_sets[i],outer_sets[j])
            inters.append(inter)
            indexes[i][j] = count
            indexes[j][i] = count
            count += 1
            intersects[i][j] = len(inter)
            intersects[j][i] = len(inter)
            jac = len(inter)/(len(outer_sets[i])+len(outer_sets[j])-len(inter)+0.0)
            print(str.format("{0}: {1}", "Jaccard", jac))
            candice[i][j] = len(inter)/(len(outer_sets[i])+0.0)
            print(str.format("{0} to {1}: {2}", names[i], names[j], candice[i][j]))
            candice[j][i] = len(inter)/(len(outer_sets[j])+0.0)
            print(str.format("{0} to {1}: {2}", names[j], names[i], candice[j][i]))


    print("Creating venn diagrams...")
    for i in range(0,len(names)-1):
        for j in range(i+1,len(names)):
            jac = intersects[i][j]/(len(outer_sets[i])+len(outer_sets[j])-intersects[i][j]+0.0)
            v = venn2(subsets=(len(outer_sets[i])-intersects[i][j],len(outer_sets[j])-intersects[i][j],intersects[i][j]),set_labels=(names[i],names[j]))
            v.get_patch_by_id('10').set_alpha(0.5)
            v.get_patch_by_id('10').set_color('orange')
            v.get_patch_by_id('01').set_alpha(0.5)
            v.get_patch_by_id('01').set_color('orange')
            v.get_patch_by_id('11').set_alpha(0.5)
            v.get_patch_by_id('11').set_color('red')
            v.get_label_by_id('10').set_text('') 
            v.get_label_by_id('01').set_text('') 
            v.get_label_by_id('11').set_text('') 
            v.get_label_by_id('10').set_size(20)
            v.get_label_by_id('01').set_size(20)
            v.get_label_by_id('11').set_size(20)
            plt.annotate('{0:,}'.format(len(outer_sets[i])-intersects[i][j]), xy = v.get_label_by_id('10').get_position(), xytext = (-30,-70), size = 'xx-large',
                ha = 'center', textcoords = 'offset points', bbox = dict(boxstyle = 'round, pad = 0.5', fc = 'lime', alpha = 0.3),
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3, rad = 0.5', color = 'gray'))
            plt.annotate('{0:,}'.format(len(outer_sets[j])-intersects[i][j]), xy = v.get_label_by_id('01').get_position(), xytext = (-30,-70), size = 'xx-large',
                ha = 'center', textcoords = 'offset points', bbox = dict(boxstyle = 'round, pad = 0.5', fc = 'lime', alpha = 0.3),
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3, rad = 0.5', color = 'gray'))
            plt.annotate('{0:,}'.format(intersects[i][j]), xy = v.get_label_by_id('11').get_position(), xytext = (-30,70), size = 'xx-large',
                ha = 'center', textcoords = 'offset points', bbox = dict(boxstyle = 'round, pad = 0.5', fc = 'lime', alpha = 0.3),
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3, rad = 0.5', color = 'gray'))
            plt.title("Venn Diagram of intersecting GIs from %s with %s (Jaccard = %.3f)" % (names[i],names[j],jac))
            plt.savefig("venn/venn_%s_%s"%(names[i],names[j]))
            plt.clf()
        # plt.show()
    print("Creating alternate jaccard bar charts...")
    for i in range(0,len(names)):
        plt.bar(range(1,len(names)),candice[i][0:i]+candice[i][(i+1):len(names)],align="center",color='orange',alpha=0.5)
        plt.title("%s's alternate Jaccard scores"%names[i])
        plt.xticks(range(1,len(names)),names[0:i]+names[(i+1):len(names)])
        plt.ylim([0,1])
        plt.savefig("alt_jaccard/alt_jaccard_%s"%names[i])
        plt.clf()
            # plt.show()

    inters3 = []
    indexes3 = [[[-1 for x in range(len(names))] for y in range(len(names))] for z in range(len(names))]
    count = 0
    for i in range(0,len(names)-2):
        for j in range(i+1,len(names)-1):
            for k in range(j+1,len(names)):
                print("merging %s %s %s..." % (names[i], names[j], names[k]))
                inters3.append(pd.merge(inters[indexes[i][j]],outer_sets[k]))
                indexes3[i][j][k] = count
                indexes3[j][i][k] = count
                print("size: {0:,}".format(len(inters3[count])))
                count += 1

    for i in tqdm(range(0,len(names)-2)):
        for j in range(i+1,len(names)-1):
            for k in range(j+1,len(names)):
                sb = (len(outer_sets[i])-intersects[i][j]-intersects[i][k]+len(inters3[indexes3[i][j][k]]),
                    len(outer_sets[j])-intersects[i][j]-intersects[j][k]+len(inters3[indexes3[i][j][k]]),
                    intersects[i][j]-len(inters3[indexes3[i][j][k]]),
                    len(outer_sets[k])-intersects[i][k]-intersects[j][k]+len(inters3[indexes3[i][j][k]]),
                    intersects[i][k]-len(inters3[indexes3[i][j][k]]),
                    intersects[j][k]-len(inters3[indexes3[i][j][k]]),
                    len(inters3[indexes3[i][j][k]]))
                v = venn3(subsets=sb,set_labels=(names[i],names[j],names[k]))
                v.get_label_by_id('100').set_text('')
                v.get_label_by_id('010').set_text('')
                v.get_label_by_id('110').set_text('')
                v.get_label_by_id('001').set_text('')
                v.get_label_by_id('101').set_text('')
                v.get_label_by_id('011').set_text('')
                v.get_label_by_id('111').set_text('')
                v.get_label_by_id('100').set_size(10)
                v.get_label_by_id('010').set_size(10)
                v.get_label_by_id('110').set_size(10)
                v.get_label_by_id('001').set_size(10)
                v.get_label_by_id('101').set_size(10)
                v.get_label_by_id('011').set_size(10)
                v.get_label_by_id('111').set_size(10)
                plt.annotate('{0:,}'.format(sb[0]), xy = v.get_label_by_id('100').get_position(), xytext = (-40,20), size = 'medium',
                ha = 'center', textcoords = 'offset points', 
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3, rad = 0.5', color = 'gray'))
                plt.annotate('{0:,}'.format(sb[1]), xy = v.get_label_by_id('010').get_position(), xytext = (50,-20), size = 'medium',
                ha = 'center', textcoords = 'offset points', 
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3, rad = 0.5', color = 'gray'))
                plt.annotate('{0:,}'.format(sb[2]), xy = v.get_label_by_id('110').get_position(), xytext = (-30,70), size = 'medium',
                ha = 'center', textcoords = 'offset points', 
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3, rad = 0.5', color = 'gray'))
                plt.annotate('{0:,}'.format(sb[3]), xy = v.get_label_by_id('001').get_position(), xytext = (40,-50), size = 'medium',
                ha = 'center', textcoords = 'offset points', 
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3, rad = 0.5', color = 'gray'))
                plt.annotate('{0:,}'.format(sb[4]), xy = v.get_label_by_id('101').get_position(), xytext = (-70,10), size = 'medium',
                ha = 'center', textcoords = 'offset points', 
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3, rad = 0.5', color = 'gray'))
                plt.annotate('{0:,}'.format(sb[5]), xy = v.get_label_by_id('011').get_position(), xytext = (70,-20), size = 'medium',
                ha = 'center', textcoords = 'offset points',
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3, rad = 0.5', color = 'gray'))
                plt.annotate('{0:,}'.format(sb[6]), xy = v.get_label_by_id('111').get_position(), xytext = (-50,20), size = 'medium',
                ha = 'center', textcoords = 'offset points',
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3, rad = 0.5', color = 'gray'))
                plt.savefig("venn/venn3_%s_%s_%s" % (names[i],names[j],names[k]))
                plt.clf()
                # plt.show()

            
    # keggseed = pd.merge(kegg, seed, how='inner', on='gi')
    # keggnog = pd.merge(kegg, eggnog, how='inner', on='gi')
    # keggcard = pd.merge(kegg, card, how='inner', on='gi')
    # kegginterpro = pd.merge(kegg, interpro, how='inner', on='gi')
    # seedeggnog = pd.merge(seed, eggnog, how='inner', on='gi')
    # seedcard = pd.merge(seed, card, how='inner', on='gi')
    # seedinterpro = pd.merge(seed, interpro, how='inner', on='gi')
    # eggnogcard = pd.merge(eggnog, card, how='inner', on='gi')
    # eggnoginterpro = pd.merge(eggnog, interpro, how='inner', on='gi')
    # cardinterpro = pd.merge(card, interpro, how='inner', on='gi')

    

    