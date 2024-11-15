#!/usr/bin/env python3
##
# @file         statsGenerator.py
#
# @author       Jakub Chlebik \n
#               Faculty of Information Technology \n
#               Brno University of Technology \n
#               xchleb07@stud.fit.vutbr.com
#
# @brief        The script file for generating statistics from the optimization results \n
#               Requires matplotlib and numpy external packages
#               [-h] [--help] for help.

# @version      0.5
#
# @date         2019-10-28 (created) \n
#               2020-04-27 (revised)
#

# IMPORTS ######################################################################
from subprocess import Popen, PIPE
import sys
import os
from datetime import datetime
import argparse
import json
import random
import statistics
from operator import itemgetter, mul

import matplotlib as mpl
mpl.use('pdf')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt
import numpy as np
# END OF IMPORTS ###############################################################

# METHODS DEFINITIONS ##########################################################

##
# @brief: Creates a JSON file object from a dictionary
#
# @param [in]       dictionary      - text output from stdout of the optimizer
# @return JSONified data            - JSON file object from a given dictionary
def createJsonFromDictionary(dictionary):
    return json.dumps(dictionary, indent=4)
# END OF createJsonFromDictionary ##############################################

##
# @brief: Creates a python dictionary from the output of the optimization algorithm 
# 
# @param [in]       optimizerOutput     - text output from stdout of the optimizer
# @return JSONified data                - python dictionary object (["Launch"]["Run"]["Result"])
def createDictFromOptimzerOutput(optimizerOutput):
    sections = zip(["Launch", "Run", "Result"], optimizerOutput.split(sep="@@@"))
    result = dict()
    for section in sections:
        result[section[0]] = dict()
        sectionLines = section[1].strip().split(sep="\n")
        for sectionPart in sectionLines:
            for keyVal in sectionPart.split(sep="$"):
                if keyVal != '' and keyVal != '\n':
                    keyVal = keyVal.split(sep=":")
                    if section[0] == "Run":
                        if not keyVal[0] in result[section[0]]:
                            result[section[0]][keyVal[0]] = list()
                        split = keyVal[1].split(sep=",")
                        if len(split) > 1:
                            result[section[0]][keyVal[0]].append([abs(float(y)) for y in split])
                        else:
                            result[section[0]][keyVal[0]].append(abs(float(keyVal[1])))
                    elif section[0] == "Result":
                        split = keyVal[1].split(sep=",")
                        if len(split) > 1:
                            result[section[0]][keyVal[0]] = [abs(float(y)) for y in split]
                        else:
                            result[section[0]][keyVal[0]] = abs(float(keyVal[1]))
                    else:
                        result[section[0]][keyVal[0]] = keyVal[1]
    if not result["Run"]:
        result.pop("Run", None)
    if not result["Launch"] or not result["Result"]:
        raise KeyError("Launch or Result missing in optimizer output") 
    return result
# END OF createDictFromOptimzerOutput ############################################

##
# @brief Creates and save the full output of the optimizer to a json file on specified path
# 
# @param [in]   outputPath          - path to save the json file to
# @param [in]   dataDictionary      - dictionary with data from the optimizer
# @param [in]   run                 - the integer run number of the optimizer in this batch
# @return Nothing
def createFullJsonOutput(outputPath, dataDictionary, run):
    outputPathPrefix = os.path.join(outputPath, str(run))
    os.makedirs(outputPathPrefix, exist_ok=True)
    jsonOutput = createJsonFromDictionary(dataDictionary)
    with open("{}/fullData.json".format(outputPathPrefix), "w+") as tempFile:
        tempFile.write(jsonOutput)
    #return jsonOutput
# END OF createFullJsonOutput ################################################################

##
# @brief Converts stats of the currently ran optimization to 4 graphs. 
# 
# @param [in]   outputPath          - path to save the plot file
# @param [in]   dataDictionary      - dictionary with data from the optimizer
# @param [in]   run                 - the run number of the optimizer in this batch
# @return Nothing
def createPlotOutput(outputPath, dataDictionary, run, outputName = "runPlot"):
    outputPathPrefix = os.path.join(outputPath, str(run))
    os.makedirs(outputPathPrefix, exist_ok=True)

    quarterOfSize = int(len(dataDictionary["gen"])*0.25)
    for j in range(0, 4):
        x = []
        data = dict() 
        data["worst"] = list()
        data["q1"] = list()
        data["median"] = list()
        data["q3"] = list()
        data["best"] = list()
        for gen in range(j*quarterOfSize, j*quarterOfSize+quarterOfSize):
            i = int(gen)
            x.append(i)
            data["q1"].append(dataDictionary["q1"][i])
            data["median"].append(dataDictionary["median"][i])
            data["q3"].append(dataDictionary["q3"][i])
            
            data["worst"].append(dataDictionary["worst"][i])
            data["best"].append(dataDictionary["best"][i])

        plt.plot(x, data["worst"], label='nejhorší')
        plt.plot(x, data["q1"], label='q1')
        plt.plot(x, data["median"], label='médián')
        plt.plot(x, data["q3"], label='q3')
        plt.plot(x, data["best"], label='nejlepší')
        plt.grid()
        plt.legend()
        plt.xlabel("Generace")
        plt.ylabel("Fitness")
        plt.savefig("{}/{}_{}.pdf".format(outputPathPrefix, outputName, j), format="pdf", bbox_inches='tight')
        plt.close()

    #x = list(range(1, len(dataDictionary["fitnessEvaluations"])+1))
    plt.plot(dataDictionary["fitnessEvaluations"], dataDictionary["best"])
    #plt.ylim(dataDictionary["fitnessEvaluations"][0], dataDictionary["fitnessEvaluations"][len(dataDictionary["fitnessEvaluations"])-1])
    plt.xlabel("Fitness evaluace")
    plt.ylabel("Nejlepší řešení")
    plt.grid()
    plt.savefig("{}/fitnessEvals.pdf".format(outputPathPrefix), format="pdf", bbox_inches='tight')
    plt.close()
# END OF createPlotOutput ################################################################

##
# @brief Create a plot graph showing the progression of the single variable during a run.\n
#        Useful for non-population based optimizers
# 
# @param [in]   outputPath          - path to save the plot file
# @param [in]   dataDictionary      - dictionary with data from the optimizer
# @param [in]   run                 - the run number of the optimizer in this batch
# @return Nothing
def createPlotOutputSingle(outputPath, dataDictionary, run, dictKeyToPlot, outputName = "runPlotBests"):
    outputPathPrefix = os.path.join(outputPath, str(run))
    os.makedirs(outputPathPrefix, exist_ok=True)

    quarterOfSize = int(len(dataDictionary["gen"])*0.25)
    for j in range(0, 4):
        x = []
        data = dict()
        data[dictKeyToPlot] = list()
        for gen in range(j*quarterOfSize, j*quarterOfSize+quarterOfSize):
            i = int(gen)
            x.append(i)
            
            data[dictKeyToPlot].append(dataDictionary[dictKeyToPlot][i])

        plt.plot(x, data[dictKeyToPlot], label='nejlepší')
        plt.grid()
        plt.legend()
        plt.xlabel("Iterace")
        plt.ylabel("Fitness")
        plt.savefig("{}/{}_{}.pdf".format(outputPathPrefix, outputName, j), format="pdf", bbox_inches='tight')
        plt.close()
# END OF createPlotOutputSingle ################################################################


##
# @brief Create an evolution progress plot, showing boxplot for each generation from median data \n
#        of all runs
# 
# @param [in]   outputPath          - path to save the plot file to
# @param [in]   dataDictionary      - dictionary with aggregated data from the optimizers (medians from all runs)
# @param [in]   outputName          - name of the plot file
# @return Nothing
def createPlotEvolutionProgress(outputPath, dataDictionary, outputName = "evolutionProgressPlot"):
    os.makedirs(outputPath, exist_ok=True)
    quarterOfSize = int(len(dataDictionary["gen"])*0.25)

    for j in range(0, 4):
        stats = []
        for gen in range(j*quarterOfSize, j*quarterOfSize+quarterOfSize):
            if gen % 2 == 1:
                i = gen + 1
            else:
                i = ""
            stats.append(
                {
                    'whislo' : dataDictionary["q1"][gen] - 1.5 * abs(dataDictionary["q3"][gen] - dataDictionary["q1"][gen]),
                    'whishi' : dataDictionary["q3"][gen] + 1.5 * abs(dataDictionary["q3"][gen] - dataDictionary["q1"][gen]),
                    'med' : dataDictionary["median"][gen],
                    'q1' : dataDictionary["q1"][gen],
                    'q3' : dataDictionary["q3"][gen],
                    'fliers' : [dataDictionary["best"][gen], dataDictionary["worst"][gen]],
                    'label' : i
                }
            )
            if stats[-1]['whislo'] < 0:
                stats[-1]['whislo'] = 0

        fig, ax1 = plt.subplots()

        ax1.set_title("Průběh jednotlivých generací {}".format(j+1))
        ax1.grid(True)
        ax1.set_xlabel("Generace")
        ax1.set_ylabel("Fitness")

        ax1.bxp(stats, patch_artist=True)

        # labels = [""] * quarterOfSize
        # fifthOfQuarter = int(quarterOfSize/5)
        # tmp = [fifthOfQuarter * i for i in range(1, 5)]

        # for i in tmp:
        #     labels[i-1] = str(i)

        # labels[0] = str(quarterOfSize*j + 1)
        # ax1.set_xticks(tmp)
        # ax1.set_xticklabels(labels)

        fig.savefig("{}/{}_{}.pdf".format(outputPath, outputName, j), format="pdf", bbox_inches='tight')
        plt.close()
# END OF createPlotEvolutionProgress ################################################################

##
# @brief Creates a plot showing the median best solution with regards to amount of fitness evaluations \n
#        performed
#
# @param [in]   outputPath          - path to save the plot file to
# @param [in]   dataDictionary      - dictionary with aggregated data from the optimizers (medians from all runs)
# @param [in]   outputName          - name of the plot file
# @return Nothing
def createPlotBestsOnFitness(outputPath, dataDictionary, outputName = "bestToFitnessPlot"):
    os.makedirs(outputPath, exist_ok=True)
    graphsNum = 2

    quarterOfSize = int(len(dataDictionary["fitnessEvaluations"])*1/graphsNum)

    for j in range(0, graphsNum):
        fitAxis = []
        data = []
        for gen in range(j*quarterOfSize, j*quarterOfSize+quarterOfSize):
            fitAxis.append(int(dataDictionary["fitnessEvaluations"][gen]))
            data.append(dataDictionary["best"][gen])

        plt.plot(fitAxis, data, label='médián nejlepších řešení')
        plt.grid()

        plt.legend()
        plt.xlabel("Počet fitness vyhodnocení")
        plt.ylabel("Fitness")
        plt.savefig("{}/{}_{}.pdf".format(outputPath, outputName, j), format="pdf", bbox_inches='tight')
        plt.close()
#END OF createPlotBestsOnFitness #######################################################################


def createPlotBestsOnFitnessMultiple(outputPath, dictList, labels, outputName = "bestToFitnessPlotMultiple"):
    os.makedirs(outputPath, exist_ok=True)
    
    minFitEv = int(dictList[0]["fitnessEvaluations"][-1])
    for i in range(0, len(dictList)):        
        if minFitEv > int(dictList[i]["fitnessEvaluations"][-1]):
            minFitEv = int(dictList[i]["fitnessEvaluations"][-1])

    for i in range(0, len(dictList)):
        dataDictionary = dictList[i]
    
        fitAxis = []
        data = []
        for gen in range(0, len(dataDictionary["fitnessEvaluations"])):
            # if  int(dataDictionary["fitnessEvaluations"][gen]) > minFitEv:
            #     fitAxis.append(minFitEv)
            #     data.append(data[-1])
            #     break
            # else:
            #     fitAxis.append(int(dataDictionary["fitnessEvaluations"][gen]))
            #     data.append(dataDictionary["best"][gen])
            fitAxis.append(int(dataDictionary["fitnessEvaluations"][gen]))
            data.append(dataDictionary["best"][gen])
        plt.plot(fitAxis, data, label=labels[i])
    plt.grid()
    plt.locator_params(axis='x', nbins=8)

    plt.legend()
    plt.xlabel("Number of fitness evaluations")
    plt.ylabel("Fitness")

    plt.savefig("{}/{}_{}.pdf".format(outputPath, outputName, 1), format="pdf", bbox_inches='tight')
    plt.close()



##
# @brief Converts the resulting stats of the whole batch into boxplots graphs
# 
# @param [in]   outputPath          - path to save the plot file
# @param [in]   results             - dictionary with fitness values of every individual in last generations in all ran optimizers in batch
# @return Nothing
def createRunsBoxplotOutput(outputPath, results, outputName = "boxplots"):
    os.makedirs(outputPath, exist_ok=True)
    plt.boxplot(results, showfliers=False, patch_artist=True, notch=True, manage_ticks=True, meanline=True)
    plt.xlabel("Běhy algoritmu")
    plt.ylabel("Fitness")

    labels = [""] * (len(results))
    tmp = [5*i for i in range(1, int(len(results)/5)+1)]

    for i in tmp:
        labels[i-1] = str(i)

    labels[0] = "1"
    plt.grid()
    plt.xticks(range(1, len(results)+1), labels)
    plt.savefig("{}/{}.pdf".format(outputPath, outputName), format='pdf', bbox_inches='tight')
    plt.close()
# END OF createRunsBoxplotOutput ################################################################

##
# @brief Creates a boxplot graph of all found solutions across the entire batch
# 
# @param [in]   outputPath          - path to save the plot file
# @param [in]   bestSolutions       - list of found solutions from last generations of all ran optimizers in batch
# @param [in]   showOutliers        - should the outliers be shown in the boxplot graph
# @return Nothing
def createRunsBestOutputsPlot(outputPath, bestSolutions, showOutliers, outputName = "bestsBoxplot", ylabel = "Fitness", xlabel = "", ticksNames = [], useNotches = False):
    os.makedirs(outputPath, exist_ok=True)

    plt.boxplot(bestSolutions, showfliers=showOutliers, patch_artist=True, notch=useNotches)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.grid(axis='y',linestyle='dashed')
    ticks, _ = plt.xticks()
    plt.xticks(ticks=ticks ,labels=ticksNames)
    plt.savefig("{}/{}.pdf".format(outputPath, outputName), format='pdf', bbox_inches='tight')
    plt.close()
# END OF createRunsBestOutputsPlot ################################################################

##
# @brief Creates a boxplot graph of fitness evaluations done across the entire batch
# 
# @param [in]   outputPath          - path to save the plot file
# @param [in]   fitnessCallsList    - list of resulting amount of fitness calls from all ran optimizers in batch
# @return Nothing
def createRunsFitnessCallsBoxplot(outputPath, fitnessCallsList):
    os.makedirs(outputPath, exist_ok=True)

    plt.boxplot(fitnessCallsList, showfliers=False, patch_artist=True)
    plt.ylabel("Počet evaluací fitness funkce")
    plt.xlabel("")
    plt.xticks(ticks=[])
    plt.savefig("{}/fitnessEvaluations.pdf".format(outputPath), format='pdf', bbox_inches='tight')
    plt.close()
# END OF createRunsFitnessCallsBoxplot ################################################################

##
# @brief Creates a folder where all statistics will be saved
# 
# @param [in]   outputPath                - a path where the statistics will be stored
# @return String                          - path to the newly created folder
def createOutputFolder(outputPath):
    now = datetime.now().replace(microsecond=0)
    outputPathPrefix = os.path.join(outputPath, now.strftime("%Y%m%d%H%M%S"))
    os.makedirs(outputPathPrefix, exist_ok=True)
    return outputPathPrefix
# END OF createRunsBestOutputsPlot ################################################################

##
# @brief Flattens a given list using recursive method
# 
# @param [in]           nestedList           - a list to be flattened
# @param [in, out]      outList              - the resulting flattened list
# @return Nothing                            - the flattened list is returned through arguments
def flatten_rec(nestedList, outList):
    for element in nestedList:
        if isinstance(element, list):
            flatten_rec(element, outList)
        else:
            outList.append(element)
# END OF flatten_rec ################################################################


def findBest(inputPath, files):
    bestFile = 0
    bestVal = 9999999999
    for i in range(0, len(files)):
        optimizerOutput = open(os.path.join(inputPath, files[i]), "r").read()
        resSection = optimizerOutput.split(sep="@@@")[-1]
        sectionLines = resSection.strip().split(sep="\n")
        for sectionLine in sectionLines:
            for keyVal in sectionLine.split(sep="$"):
                keyVal = keyVal.split(sep=":")
                if keyVal[0] == "fitness":
                    val = float(keyVal[1])
                    if val < bestVal:
                        bestVal = val
                        bestFile = i
                        break
            break
    return bestFile

##
# @brief Core of the script. Creates statistical graphs from the given folder and outputs the \n
#        results.
# 
# @param [in] args           - arguments given to the script
# @return Nothing            - the flattened list is returned through arguments
def generateStatistics(args, createGraphs, files = []):
    if createGraphs:
        outputPath = createOutputFolder(args.outputFolder)

    everyRunStats = True
    isEvo = bool(args.isPop)

    lastGenFits = []
    results = []
    entireRunFits = []
    
    timeNeeded = []

    bests = []
    medians = []
    q1s = []
    q3s = []
    worsts = []
    fitnessEvals = []

    if not files:
        files = [f for f in next(os.walk(args.inputFolder))[2] if f.endswith(".optOut")]

    longestRunSize = 0
    i = 0
    for f in files:
        optimizerOutput = open(os.path.join(args.inputFolder, f), "r").read()
        dataDictionary = createDictFromOptimzerOutput(optimizerOutput)
        if createGraphs:
            createFullJsonOutput(outputPath, dataDictionary, i)
        if "Run" in dataDictionary :

            if everyRunStats and createGraphs:
                createPlotOutput(outputPath, dataDictionary["Run"], i, "runPlot")
                createPlotOutputSingle(outputPath, dataDictionary["Run"], i, "best", "runPlotBestsOnly")
    
            lastPopIndex = len(dataDictionary["Run"]["populationDump"])-1

            #boxplots of the last generations
            lastGenFits.append(dataDictionary["Run"]["populationDump"][lastPopIndex])

            #boxplot of best results
            results.append(dataDictionary["Result"]["fitness"])

            #boxplots of the populations for the entire run
            flatList = []
            flatten_rec(dataDictionary["Run"]["populationDump"], flatList)
            entireRunFits.append(flatList)

            #boxplot of time needed for evaluations
            timeNeeded.append(dataDictionary["Result"]["time[uS]"]/1000000)
            
            #data collection for the evolution progress boxplots
            bests.append(dataDictionary["Run"]["best"])
            worsts.append(dataDictionary["Run"]["worst"])
            medians.append(dataDictionary["Run"]["median"])
            q1s.append(dataDictionary["Run"]["q1"])
            q3s.append(dataDictionary["Run"]["q3"])
            fitnessEvals.append(dataDictionary["Run"]["fitnessEvaluations"])

            if len(dataDictionary["Run"]["best"]) > longestRunSize:
                longestRunSize = len(dataDictionary["Run"]["best"])
                longestRunNum = i
        else:
            print("Verbose flag was not set for optimizers, cannot generate statistics from file " + f, file=sys.stderr)
        i += 1


    aggregatedData = {}
    aggregatedData["gen"] = []
    aggregatedData["fitnessEvaluations"] = []
    aggregatedData["best"] = []
    aggregatedData["median"] = []
    aggregatedData["q1"] = []
    aggregatedData["q3"] = []
    aggregatedData["worst"] = []

    tempListVals = []
    tempListVals.append([])
    tempListVals.append([])
    tempListVals.append([])
    tempListVals.append([])
    tempListVals.append([])

    for gen in range(0, longestRunSize):
        for run in range(0, len(bests)):
            if len(bests[run])-1 < gen: # run shorter than max, fill rest of values with the last
                safeIndex = len(bests[run])-1
            else:
                safeIndex = gen

            tempListVals[0].append(bests[run][safeIndex])
            tempListVals[1].append(medians[run][safeIndex])
            tempListVals[2].append(q1s[run][safeIndex])
            tempListVals[3].append(q3s[run][safeIndex])
            tempListVals[4].append(worsts[run][safeIndex])

        aggregatedData["gen"].append(gen)
        aggregatedData["best"].append(statistics.median(tempListVals[0]))
        aggregatedData["median"].append(statistics.median(tempListVals[1]))
        aggregatedData["q1"].append(statistics.median(tempListVals[2]))
        aggregatedData["q3"].append(statistics.median(tempListVals[3]))
        aggregatedData["worst"].append(statistics.median(tempListVals[4]))
        aggregatedData["fitnessEvaluations"].append(fitnessEvals[longestRunNum][gen])

        tempListVals[0] = []
        tempListVals[1] = []
        tempListVals[2] = []
        tempListVals[3] = []
        tempListVals[4] = []
    if createGraphs:
        if isEvo:
            createPlotEvolutionProgress(outputPath,  aggregatedData, "evolutionProgress")
        createPlotBestsOnFitness(outputPath,     aggregatedData, "bestsToFitness")

        createRunsBoxplotOutput(outputPath, lastGenFits,    "lastGenBoxplots")
        createRunsBoxplotOutput(outputPath, entireRunFits,  "entirePopulationBoxplot")

        createRunsBestOutputsPlot(outputPath,   results,    False,  "bestsBoxplot")
        createRunsBestOutputsPlot(outputPath,   results,    True,   "bestsBoxplot_WithOutliers")
        
        createRunsBestOutputsPlot(outputPath,   timeNeeded,    True,   "timeBoxplot_WithOutliers", "Čas [s]")
    return (results, timeNeeded, list(map(itemgetter(-1), fitnessEvals)),  aggregatedData)
# END OF generateStatistics ######################################################

# END OF METHOD DEFINITIONS ######################################################



##
# @brief: Main entry point of the script. Takes the path to folder containing the results of \n
#         optimizers and the path to a folder where to store the statistics. This script \n
#         creates statistical graphs from the results of the optimizers.
# 
def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    #group2 = parser.add_argument_group('comp')

    group.add_argument('-i', '--inputFolder', help="Path to folder containing the outputs of optimization runs. Every file with .optOut extension is used.")
    group.add_argument('-m', '--multipleFolders', help="Many different paths to folders containing statistics to aggregate.", nargs='+')
    group.add_argument('-c', '--compareFolderResults', help="TODO.", nargs='+')

    parser.add_argument('-o', '--outputFolder', help="Path to a folder where to store the statistical results", required=True)
    parser.add_argument('-p', '--isPop', help="Are the results from population based algorithm", type=int, required=True)
    
    parser.add_argument('-n', '--names', help="TODO", nargs='+')
    

    args = parser.parse_args()
    try:
        if args.compareFolderResults:
            resList = []
            timesList = []
            methodNames = []
            fitnessEvals = []
            tmp = []
            for i in range(0, len(args.compareFolderResults)):
                Fs = []
                evalFactor = 1
                args.inputFolder = args.compareFolderResults[i]
                for root, _, files in os.walk(args.inputFolder, topdown=False):
                    if files:
                        if len(files) > 1:
                            evalFactor = len(files)
                            Fs.append(os.path.join(root.split('/')[-1], files[findBest(root, files)]))
                        else:
                            Fs.append(os.path.join(root.split('/')[-1], files[0]))

                methodName = args.names[i]
                (res, timeTaken, evalsNeeded, aggregateData) = generateStatistics(args, False, Fs)
                resList.append(res)
                tmp.append(aggregateData)
                fitnessEvals.append([fitEval * evalFactor for fitEval in evalsNeeded])
                timesList.append([time/3600 for time in timeTaken])
                methodNames.append(methodName)
            outputPath = createOutputFolder(args.outputFolder)
            #createPlotBestsOnFitnessMultiple(outputPath,     tmp,  methodNames, "bestsToFitness")
            createRunsBestOutputsPlot(outputPath, resList, False, "solutionsPlotsComparasion", "Fitness", "", methodNames, False)
            createRunsBestOutputsPlot(outputPath, timesList, False, "runTimesPlotsComparasion", "Time [h]", "", methodNames, False)
            createRunsBestOutputsPlot(outputPath, fitnessEvals, False, "fitEvalsPlotsComparasion", "Fitness Evaluations", "", methodNames, False)
            print(resList, timesList)
        elif args.multipleFolders:
            resList = []
            timesList = []
            methodNames = []
            fitnessEvals = []
            Fs = []
            for i in range(0, len(args.multipleFolders)):
                args.inputFolder = args.multipleFolders[i]                

                methodName = args.inputFolder.split('/')[-1]
                (res, timeTaken, evalsNeeded, _) = generateStatistics(args, False)
                resList.append(res)
                fitnessEvals.append(evalsNeeded)
                timesList.append(timeTaken)
                methodNames.append(methodName)
            outputPath = createOutputFolder(args.outputFolder)
            createRunsBestOutputsPlot(outputPath, resList, False, "solutionsPlotsComparasion", "Fitness", "", methodNames, False)
            createRunsBestOutputsPlot(outputPath, timesList, False, "runTimesPlotsComparasion", "Čas [s]", "", methodNames, True)
        else:
            generateStatistics(args, True)
    except Exception as err:
        print(err)


if __name__ == "__main__":
    main()
# END of MAIN ###################################################################