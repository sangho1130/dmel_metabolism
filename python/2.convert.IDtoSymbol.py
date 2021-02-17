#!/usr/bin/python

import argparse
import sys

def chg_Id2Sym(exprArg):
    annoTable = open('fbgn_annotation_ID_fb_2019_03_filtered.txt', 'r')
    tableLines = annoTable.readlines()
    annoTable.close()

    annoDict = dict()
    redunList = list(); redunDict = dict()
    count = 0
    for tableLine in tableLines[1:]:
        tableLine = tableLine.rstrip().split('\t')
        annoDict[tableLine[2]] = tableLine[0]

        if len(tableLine[3]) == 0:  continue
        secondIds = tableLine[3].strip('"').split(',')
        if len(secondIds) >= 1:
            for secondId in secondIds:
                if secondId in redunList:
                    redunDict[secondId].append(tableLine[0])
                    continue
                elif annoDict.has_key(secondId):
                    redunList.append(secondId)
                    if redunDict.has_key(secondId):
                        redunDict[secondId].append(tableLine[0])
                        redunDict[secondId].append(annoDict[secondId])
                    else:
                        redunDict[secondId] = [tableLine[0], annoDict[secondId]]
                    annoDict.pop(secondId, None)
                    continue
                else:
                    annoDict[secondId] = tableLine[0]
    
    openExpr = open(exprArg, 'r')
    exprLines = openExpr.readlines()
    openExpr.close()
    
    outwrite = open(exprArg[:-3] + 'symbol.txt', 'w')
    outwrite.write(exprLines[0])
    writeDict = dict()
    for exprLine in exprLines[1:]:
        exprLine = exprLine.rstrip().split('\t')
        if exprLine[0] in ['EGFP', 'gal4']:
            continue
        elif not annoDict.has_key(exprLine[0]):
            if exprLine[0] in redunList:
                print 'in redunList', exprLine[0]
            else:
                print 'missing', exprLine[0]
            continue
        elif writeDict.has_key(annoDict[exprLine[0]]):
            print 'redun', exprLine[0], annoDict[exprLine[0]], ';'
            writeDict.pop(annoDict[exprLine[0]], None)
            continue

        writeLine = annoDict[exprLine[0]] + '\t' + '\t'.join(exprLine[1:]) + '\n'
        writeDict[annoDict[exprLine[0]]] = writeLine
    print len(writeDict.keys())
    for akey in sorted(writeDict.keys()):
        outwrite.write(writeDict[akey])

    outwrite.close()


def main(args):
    chg_Id2Sym(args.Expr)

if __name__ == '__main__':
    parser = argparse.ArgumentParser('')
    parser.add_argument('-e', '--Expr', help = '', required = True)
    args = parser.parse_args()
    main(args)
