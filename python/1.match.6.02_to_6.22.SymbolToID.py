#!/home/sangho/miniconda3/envs/py2/bin/python

#
#

import argparse

def chg_Sym2Id(exprArg, outputArg, gtfArg):
    openGtf = open(gtfArg, 'r')
    gtfLines = openGtf.readlines()
    openGtf.close()
    gtfDict = dict()
    for gtfLine in gtfLines:
        if gtfLine.rstrip().split()[2] == 'gene':
            gtfLine = gtfLine.rstrip().split('\t')
            sym = gtfLine[-1].split()[3].strip('";')
            id = gtfLine[-1].split()[1].strip('";')
            gtfDict[sym] = id

    outwrite = open(outputArg, 'w')
    
    openExpr = open(exprArg, 'r')
    exprLines = openExpr.readlines()
    openExpr.close()
    outwrite.write(exprLines[0])

    writeDict = dict()
    for exprLine in exprLines[1:]:
        writeLine = gtfDict[exprLine.rstrip().split()[0]] + '\t' + '\t'.join(exprLine.rstrip().split()[1:]) + '\n'
        writeDict[gtfDict[exprLine.rstrip().split()[0]]] = writeLine
    
    for akey in sorted(writeDict.keys()):
        outwrite.write(writeDict[akey])
    outwrite.close()


def main(args):
    chg_Sym2Id(args.Expr, args.Output, args.Gtf)

if __name__ == '__main__':
    parser = argparse.ArgumentParser('')
    parser.add_argument('-e', '--Expr', help = 'Expression table with fist coloumn as symbols', required = True)
    parser.add_argument('-o', '--Output', help = 'Output file path and name', required = True)
    parser.add_argument('-g', '--Gtf', help = 'GTF annotation for expression table', required = False,
                        default = 'Drosophila_melanogaster.BDGP6.withLincrnas.gtf')
    args = parser.parse_args()
    main(args)
