#!/usr/bin/python3

import boto3
import datetime
import os
import sys
import math
import copy
import itertools
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
    total = sum(readLengths)
    return(total, percentiles)


def getReferenceStats(referenceFilePath):
    referenceLengths = []
    with open(referenceFilePath) as fp:
        for line in fp:
            if line[0] == '>':
                continue
            referenceLengths.append(len(line))
    
    total = sum(referenceLengths)
    return(total)

def processQuastReport(reportSummaries, overrideParams, quastReportFilePath):
    keysToExtract = {
        "Genome fraction (%)",
        "NG50",
        "# misassemblies",
        "# mismatches per 100 kbp",
        "# indels per 100 kbp"
    }
    configParams = {}
    for p in overrideParams:
        (k, v) = p.split()
        k = k.strip('-')
        configParams[k] = v

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


def getParamLists(paramDict):
    """
    Values in paramDict are arrays. This method explodes and flattens as follows ...
    E.g. paramDict = {'x': [1,2,3], 'y': [4,5]}
    is converted to
    [['--x 1', '--x 2', '--x 3'], ['--y 4', '--y 5']]
    """
    lists = []
    for (k, vArr) in paramDict.items():
        values = []
        for v in vArr:
            values.append("--{} {}".format(k, v))
        lists.append(values)
    return lists


def runThisConfig(shastaBinaryPath,
    inputFilePath,
    referenceFilePath,
    workingDirPath,
    overrideParams,
    reportSummaries,
    totalAssemblyTime,
    totalQATime):

    assemblyDirPath = '{}/assembly'.format(workingDirPath)
    assemblyPath = '{}/Assembly.fasta'.format(assemblyDirPath)

    os.system('rm -rf {}'.format(assemblyDirPath))

    configFlags = ' '.join(overrideParams)
    command = """
    {} --input {} --Align.alignMethod 3 --assemblyDirectory {} {}
    """.format(
        shastaBinaryPath,
        inputFilePath,
        assemblyDirPath,
        configFlags
    )

    begin = datetime.datetime.now()
    os.system(command)
    end = datetime.datetime.now()
    totalAssemblyTime = totalAssemblyTime + (end - begin).seconds

    quastDirPath = '{}/quast_results'.format(workingDirPath)
    quast_command = """
    quast.py --large --threads {} --min-identity 80 -r {} -o {} {}
    """.format(
        multiprocessing.cpu_count(),
        referenceFilePath,
        quastDirPath,
        assemblyPath
    )
    begin = datetime.datetime.now()
    os.system(quast_command)
    end = datetime.datetime.now()
    totalQATime = totalQATime + (end - begin).seconds

    quastReportFilePath = '{}/quast_results/report.tsv'.format(workingDirPath)
    processQuastReport(reportSummaries, overrideParams, quastReportFilePath)


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

    (readTotal, percentiles) = getReadStats(inputFilePath)
    referenceTotal = getReferenceStats(referenceFilePath)

    coverage = float(readTotal) / float(referenceTotal)
    # Use coverage for coming up with parameter ranges for some params!

    paramUniverse = {}
    paramUniverse['Reads.minReadLength'] = [
        math.floor(percentiles[0]), # 10th percentile
        math.floor(percentiles[1]), # 20th percentile
        math.floor(percentiles[2]), # 30th percentile
        10000, # Default value
    ]
    paramUniverse['Kmers.k'] = [
        9,
        10, # Default value
    ]

    # This can be a function of read lengths & Kmers.k??
    paramUniverse['MinHash.m'] = [
        3,
        4, # Default value
        5,
    ]
    paramUniverse['MinHash.alignmentCandidatesPerRead'] = [
        15,
        20, # Default value
        25,
    ]
    
    # Generate all possible combinations of param values.
    overrideParamCombos = itertools.product(*getParamLists(paramUniverse))
    
    totalAssemblyTime = 0
    totalQATime = 0
    for overrideParams in overrideParamCombos:
        runThisConfig(
            shastaBinaryPath,
            inputFilePath,
            referenceFilePath,
            'workingDir',
            overrideParams,
            reportSummaries,
            totalAssemblyTime,
            totalQATime
        )

    reportFilePath = '{}-summary.tsv'.format(os.path.basename(inputFilePath).split('.')[0])            
    writeSummaries(reportSummaries, reportFilePath)

    print("Total Assembly time = {} s".format(totalAssemblyTime))
    print("Total QA time = {} s".format(totalQATime))

    # Upload summary to s3.
    # s3_client = boto3.client('s3')
    # try:
    #     response = s3_client.upload_file(reportFilePath, 'czi.bhal-public', os.path.basename(reportFilePath))
    # except ClientError as e:
    #     print(e)
    
    return


if __name__ == '__main__':
    main()