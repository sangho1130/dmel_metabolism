#!/usr/bin/python2

# 2018-05-01; version 1.0.0
# Merge drop-seq libraries; set genes to the smallest library

import argparse
import sys

def mergeLibs(exprListArg, outputArg, wholeFlag, annoArg):
    if annoArg == 1:
        header = 'Symbol'
    elif annoArg == 2:
        header = 'ID\tSymbol'
    exprDict = dict()
    sampleCountDict = dict()
    for exprtable in exprListArg:
        openExpr = open(exprtable, 'r')
        exprLines = openExpr.readlines()
        openExpr.close()
        
        if annoArg == 1:
            header += '\t' + '\t'.join(exprLines[0].rstrip().split()[1:])
            sampleCountDict[exprtable] = len(exprLines[0].rstrip().split()[1:])
        elif annoArg == 2:
            header += '\t' + '\t'.join(exprLines[0].rstrip().split()[2:])
            sampleCountDict[exprtable] = len(exprLines[0].rstrip().split()[2:])
        tmpDict = dict()
        for exprLine in exprLines[1:]:
            if annoArg == 1:
                tmpDict[exprLine.rstrip().split()[0]] = '\t'.join(exprLine.rstrip().split()[1:])
            elif annoArg == 2:
                tmpDict['\t'.join(exprLine.rstrip().split()[:2])] = '\t'.join(exprLine.rstrip().split()[2:])
        exprDict[exprtable] = tmpDict
    header += '\n'
    
    if wholeFlag:#use all genes
        geneList = set()
        for exprtable in exprListArg:
            geneList = geneList|set(exprDict[exprtable].keys())
        print len(geneList), 'genes detected'
    else:
        geneList = set()
        firstDict = True
        for exprtable in exprListArg:
            if firstDict:
                geneList = set(exprDict[exprtable].keys())
                firstDict = False
            else:
                geneList = geneList&set(exprDict[exprtable].keys())
        print len(geneList), 'genes detected'
    geneList = list(geneList)
    geneList.sort()

    outwrite = open(outputArg, 'w')
    outwrite.write(header)
    for gene in geneList:
        writeLine = gene
        for exprtable in exprListArg:
            if exprDict[exprtable].has_key(gene):
                writeLine += '\t' + exprDict[exprtable][gene]
            else:
                writeLine += '\t' + '\t'.join(['0']*sampleCountDict[exprtable])
        outwrite.write(writeLine + '\n')
    outwrite.close()



def main(args):
    mergeLibs(args.Expr, args.Output, args.whole, args.Anno)

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Merge drop-seq libraries; set genes to the smallest library')
    parser.add_argument('-e', '--Expr', nargs='+',  help = 'Expression tables, multiple files in ordered and space-deliminated', required = True)
    parser.add_argument('-o', '--Output', help = 'An output file path and name', required = True)
    parser.add_argument('--use-all', dest = 'whole', help = 'Set gene list to the largest gene set (default is set to the smallest gene set)', action = 'store_true')
    parser.add_argument('-a', '--Anno', help = 'number of annotations; default = 1', required = False, default = 1, type = int)
    parser.set_defaults(whole = False)
    args = parser.parse_args()
    main(args)
