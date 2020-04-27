#!/usr/bin/python3

import os
import sys
import math
import copy
import multiprocessing
import numpy as np

helpMessage = """
This script will take as input
1. A file containing reads in the fasta format
2. A file containing the reference assembly
and it will run a few experiments to find a good Shasta config.

Usage:
    SuggestConfig.py input.fasta reference.fasta
"""

def getReadStats(inputFilePath):
    readLengths = []
    
    with open(inputFilePath) as fp:
        for line in fp:
            line = line.strip()
            if line[0] == '>':
                continue
            readLengths.append(len(line))

    percentiles = np.percentile(readLengths, [10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    return(percentiles)


def getReferenceStats(referenceFilePath):
    readLengths = []
    with open(referenceFilePath) as fp:
        for line in fp:
            if line[0] == '>':
                continue
            readLengths.append(len(line))
    
    total = sum(readLengths)
    return(total)

def processQuastReport(reportSummaries, configParams, quastReportFilePath):
    keysToExtract = {
        "Genome fraction (%)",
        "NG50",
        "# misassemblies",
        "# mismatches per 100 kbp",
        "# indels per 100 kbp"
    }
    reportDict = copy.deepcopy(configParams)
    with open(quastReportFilePath) as fp:
        for line in fp:
            line = line.strip()
            (key, value) =  line.split("\t")
            if key not in keysToExtract:
                continue
            reportDict[key] = float(value)
    
    reportSummaries.append(reportDict)
    return

def writeSummaries(reportSummaries, reportFilePath):
    keys = reportSummaries[0].keys()
    keys = sorted(keys)
    header = '\t'.join(keys)
    with open(reportFilePath, 'w') as fp:
        fp.write('{}\n'.format(header))

        for reportSummary in reportSummaries:
            row = []
            for key in keys:
                row.append(reportSummary[key])
            fp.write('{}\n'.format('\t'.join(str(entry) for entry in row)))
    return


def runThisConfig(shastaBinaryPath,
    inputFilePath,
    referenceFilePath,
    workingDirPath,
    configParams,
    reportSummaries):

    assemblyDirPath = '{}/assembly'.format(workingDirPath)
    assemblyPath = '{}/Assembly.fasta'.format(assemblyDirPath)

    os.system('rm -rf {}'.format(assemblyDirPath))

    command = """
    {} --input {} --assemblyDirectory {}  --Reads.minReadLength {} --Align.minAlignedMarkerCount {}
    """.format(
        shastaBinaryPath,
        inputFilePath,
        assemblyDirPath,
        configParams['Reads.minReadLength'],
        configParams['Align.minAlignedMarkerCount']
    )
    os.system(command)

    quastDirPath = '{}/quast_results'.format(workingDirPath)
    quast_command = """
    quast.py --large --threads {} --min-identity 80 -r {} -o {} {}
    """.format(
        multiprocessing.cpu_count(),
        referenceFilePath,
        quastDirPath,
        assemblyPath
    )
    os.system(quast_command)

    quastReportFilePath = '{}/quast_results/report.tsv'.format(workingDirPath)
    processQuastReport(reportSummaries, configParams, quastReportFilePath)


def main():
    if not len(sys.argv)==4:
        print(helpMessage)
        sys.exit(1)

    shastaBinaryPath = sys.argv[1]
    inputFilePath = sys.argv[2]
    referenceFilePath = sys.argv[3]

    if not os.path.isfile(shastaBinaryPath):
        print("Shasta binary not found at {}".format(shastaBinaryPath))
        sys.exit(1)

    if not os.path.isfile(inputFilePath):
        print("File path {} does not exist. Exiting...".format(inputFilePath))
        sys.exit(1)

    if not os.path.isfile(referenceFilePath):
        print("File path {} does not exist. Exiting...".format(referenceFilePath))
    
    # List of Quast report summaries along with their corresponding shasta config overrides.
    reportSummaries = []

    readStats = getReadStats(inputFilePath)
    # referenceStats = getReferenceStats(referenceFilePath)

    paramUniverse = {}
    paramUniverse['Reads.minReadLength'] = [
        math.floor(readStats[0]), # 10th percentile
        math.floor(readStats[1]), # 20th percentile
        math.floor(readStats[2]), # 30th percentile
    ]
    
    for rl in paramUniverse['Reads.minReadLength']:
        for markerCountDivisor in [90, 100, 110]: 
            configParams = {}
            configParams['Reads.minReadLength'] = rl
            configParams['Align.minAlignedMarkerCount'] = math.floor(rl/markerCountDivisor)

            runThisConfig(
                shastaBinaryPath,
                inputFilePath,
                referenceFilePath,
                'workingDir',
                configParams,
                reportSummaries
            )

    reportFilePath = 'summary.tsv'            
    writeSummaries(reportSummaries, reportFilePath)



if __name__ == '__main__':
    main()