import matplotlib as mpl
mpl.use('pdf')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['axes.axisbelow'] = True

import matplotlib.pyplot as plt
import numpy as np


colors = ['indianred', 'khaki']

flower = True
blob = True
blobTimes = True
scaling = False

if blobTimes:
    x = ["Farmer6x6", "Island6x6", "Farmer20x6", "Island20x6"]

    timesList = [
        [1.0779471813888888, 0.6732699372222222, 8.005410630277778, 0.45707872833333335, 0.8889327719444444, 6.430171958333333, 7.8318886319444445, 0.6789095716666667, 0.7362457786111111, 8.004701922499999, 1.4506682969444444, 0.9197000916666667, 0.8233934086111111, 7.0978730047222225, 8.004653892222223], 
        [7.886087108888889, 1.4293157419444444, 1.4476386358333333, 1.933451123888889, 1.9268362322222221, 2.9692531547222223, 8.000395083888888, 1.7535114838888888, 2.5289313475000004, 6.2699973375, 3.8570941458333334, 1.849712545, 1.720579943611111, 1.4628383969444445, 1.5606562791666667],
        [5.173041295833333, 5.500502944722222, 5.028428783055555, 5.831511598333333, 4.413556479166667, 5.217050149166667, 3.4919168380555554, 3.258855526111111, 4.0401766497222225, 4.960628355277778, 5.3294922244444445, 5.411930111388889, 3.644810633055555, 4.937645193611111, 4.2754881475], 
        [8.024208969444445, 8.020690788888889, 8.015354873888889, 8.018357222222223, 8.010404246111111, 8.003652144166667, 8.006526471111112, 8.017193805833333, 8.018142248055556, 8.0076924375, 8.008868800277778, 8.010977135833333, 8.021138545000001, 8.015348444166667, 8.01581799611111]
    ]

    bplot = plt.boxplot(timesList, showfliers=False, patch_artist=True)
    plt.axvline(x=2.5, color='k', linestyle=':')
    plt.grid(axis='y',linestyle='dashed')
    plt.ylabel('Time [h]') 
    ticks, _ = plt.xticks()
    plt.xticks(ticks=ticks ,labels=x)

    bplot['boxes'][0].set(facecolor=colors[0])
    bplot['medians'][0].set(color='black')
    bplot['boxes'][1].set(facecolor=colors[1])
    bplot['medians'][1].set(color='black')
    bplot['boxes'][2].set(facecolor=colors[0])
    bplot['medians'][2].set(color='black')
    bplot['boxes'][3].set(facecolor=colors[1])
    bplot['medians'][3].set(color='black')

    plt.savefig("timesTaken.pdf", bbox_inches='tight')
    plt.close()
#FLOWER ALL
if flower:
    x = ["Farmer6x6", "Island6x6", "Farmer20x6", "Island20x6"]

    timesList = [
        [1.0779471813888888, 0.6732699372222222, 8.005410630277778, 0.45707872833333335, 0.8889327719444444, 6.430171958333333, 7.8318886319444445, 0.6789095716666667, 0.7362457786111111, 8.004701922499999, 1.4506682969444444, 0.9197000916666667, 0.8233934086111111, 7.0978730047222225, 8.004653892222223], 
        [7.886087108888889, 1.4293157419444444, 1.4476386358333333, 1.933451123888889, 1.9268362322222221, 2.9692531547222223, 8.000395083888888, 1.7535114838888888, 2.5289313475000004, 6.2699973375, 3.8570941458333334, 1.849712545, 1.720579943611111, 1.4628383969444445, 1.5606562791666667],
        [5.173041295833333, 5.500502944722222, 5.028428783055555, 5.831511598333333, 4.413556479166667, 5.217050149166667, 3.4919168380555554, 3.258855526111111, 4.0401766497222225, 4.960628355277778, 5.3294922244444445, 5.411930111388889, 3.644810633055555, 4.937645193611111, 4.2754881475], 
        [8.024208969444445, 8.020690788888889, 8.015354873888889, 8.018357222222223, 8.010404246111111, 8.003652144166667, 8.006526471111112, 8.017193805833333, 8.018142248055556, 8.0076924375, 8.008868800277778, 8.010977135833333, 8.021138545000001, 8.015348444166667, 8.01581799611111]
    ]

    resList = [
        [51.0, 500.0, 65.0, 52.0, 151.042, 61.0, 79.0, 63.0, 81.0, 41.0, 90.0, 164.0, 109.0, 105.0, 48.0], 
        [494.792, 646.875, 538.0, 576.042, 840.625, 390.833, 201.5, 213.0, 153.125, 160.0, 147.917, 108.333], 
        [958.1, 769.575, 1354.5, 1399.43, 1834.69, 1547.49, 1673.19, 1401.75, 830.85, 685.875, 1603.88, 1258.43, 1175.04, 1635.38, 1542.64], 
        [1618.2, 1550.4, 1595.65, 1697.03, 1463.75, 1440.44, 1587.16, 1588.38, 1568.25, 1488.53, 1645.88, 1302.9, 1344.25, 1800.56, 1487.5]
    ]

    bplot = plt.boxplot(resList, showfliers=False, patch_artist=True)
    plt.axvline(x=2.5, color='k', linestyle=':')
    plt.grid(axis='y',linestyle='dashed')
    plt.ylabel('Fitness') 
    ticks, _ = plt.xticks()
    plt.xticks(ticks=ticks ,labels=x)

    bplot['boxes'][0].set(facecolor=colors[0])
    bplot['medians'][0].set(color='black')
    bplot['boxes'][1].set(facecolor=colors[1])
    bplot['medians'][1].set(color='black')
    bplot['boxes'][2].set(facecolor=colors[0])
    bplot['medians'][2].set(color='black')
    bplot['boxes'][3].set(facecolor=colors[1])
    bplot['medians'][3].set(color='black')

    ax = plt.axes([.17, .60, .2, .2])
    bplot = plt.boxplot(resList[0], showfliers=False, patch_artist=True)
    plt.title('Farmer6x6')
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    plt.legend(frameon=False)
    plt.grid(axis='y',linestyle='dashed')

    bplot['boxes'][0].set(facecolor=colors[0])
    bplot['medians'][0].set(color='black')

    plt.savefig("flowerRes.pdf", bbox_inches='tight')
    plt.close()


#BLOB
if blob:
    x = ["Farmer6x6", "Island6x6"]
    data = [[0,0,8.3,0,0,10.83,8,0,0,4,0,0,0,10,8], [0,0,0,0,0,0,0,0,0,0,0,0,0, 4.5, 12.5]]
    scsRate = [9/15*100, 15/15*100]

    #plt.set_title('Basic Plot')
    bplot = plt.boxplot(data, showfliers=False, patch_artist=True)
    plt.grid(axis='y',linestyle='dashed')
    plt.ylabel('Fitness') 
    ticks, _ = plt.xticks()
    plt.xticks(ticks=ticks ,labels=x)

    bplot['boxes'][0].set(facecolor=colors[0])
    bplot['medians'][0].set(color='black')
    bplot['boxes'][1].set(facecolor=colors[1])
    bplot['medians'][1].set(color='black')

    a = plt.axes([.55, .5, .3, .3])
    bar = plt.bar(x, scsRate, edgecolor = "k")
    plt.title('Success Rate')
    plt.rc('axes', axisbelow=True)
    #plt.xticks(ticks2, ticks2)
    plt.legend(frameon=False)
    plt.grid(axis='y',linestyle='dashed')
    #plt.ylabel('Percentage') 
    #plt.show()
    bar[0].set_color(colors[0])
    bar[1].set_color(colors[1])
    plt.savefig("6.pdf", bbox_inches='tight')

    plt.close()

    data = [ 
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,5.93,10.2,6.1875,0,3.2625,1.1625,11.375,5.875,0,7.275,2.2,10.9125,10.875,18.8]
        ] 
    scsRate = [15/15*100, 3/15*100]

    x = ["Farmer20x6", "Island20x6"]
    bplot = plt.boxplot(data, showfliers=False, patch_artist=True)
    plt.grid(axis='y',linestyle='dashed')
    plt.ylabel('Fitness') 
    ticks, _ = plt.xticks()
    plt.xticks(ticks=ticks ,labels=x)

    bplot['boxes'][0].set(facecolor=colors[0])
    bplot['medians'][0].set(color='black')
    bplot['boxes'][1].set(facecolor=colors[1])
    bplot['medians'][1].set(color='black')

    ax = plt.axes([.17, .5, .3, .3])
    bar = plt.bar(x, scsRate, edgecolor = "k")
    plt.title('Success Rate')
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    plt.legend(frameon=False)
    plt.grid(axis='y',linestyle='dashed')
    #plt.ylabel('Percentage') 
    #plt.show()
    bar[0].set_color(colors[0])

    bar[1].set_color(colors[1])
    plt.savefig("20.pdf", bbox_inches='tight')




#SCALING
if scaling:
    data = [10.2561, 6.7482, 2.79, 1.66641, 1.4049, 1.29487, 1.20264]
    ticks = [1, 2, 6, 12, 16, 20, 24]

    plt.plot(ticks, data, marker='o', linestyle='-', color=colors[0]) 
    plt.xlabel('Cores')
    plt.ylabel('Time [s]') 
    plt.xticks(ticks, ticks)
    plt.legend(frameon=False)
    plt.grid()


    data2 = [1.66641, 1.4049, 1.29487, 1.20264]
    ticks2 = [12, 16, 20, 24]

    a = plt.axes([.55, .5, .3, .3])
    plt.plot(ticks2, data2, marker='o', linestyle='-', color='b')
    plt.title('Detail')

    plt.xticks(ticks2, ticks2)
    plt.legend(frameon=False)
    plt.grid()
    plt.savefig("scalingZoom.pdf")
    plt.close()

