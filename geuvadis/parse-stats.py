#!/home/ruolin/Software/anaconda3/bin/python
import sys
import argparse
import re

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-f", nargs="+", default = [], required=True)
    args = parser.parse_args()
    files = args.f

    TEXON = re.compile(r"^.*Missed exons:\s+(\d+)\/(\d+)*.")
    FEXON = re.compile(r"^.*Novel exons:\s+(\d+)\/(\d+)*.")

    TINTRON = re.compile(r"^.*Missed introns:\s+(\d+)\/(\d+)*.")
    FINTRON = re.compile(r"^.*Novel introns:\s+(\d+)\/(\d+)*.")

    TLOCI = re.compile(r"^.*Missed loci:\s+(\d+)\/(\d+)*.")
    FLOCI = re.compile(r"^.*Novel loci:\s+(\d+)\/(\d+)*.")

    TTX = re.compile(r"^.*Matching intron chains:\s+(\d+)")
    QTX = re.compile(r"^#.*Query mRNAs\s+:\s+(\d+)*.")
    RTX = re.compile(r"^#.*Reference mRNAs\s+:\s+(\d+)*.")

    exon_sen = []
    exon_tpp = []

    intron_sen = []
    intron_tpp = []

    loci_sen = []
    loci_tpp = []

    tx_sen = []
    tx_tpp = []

    for f in files:
        with open(f, 'r') as fh:
            tx_tp = 0.0
            tx_nref = 0.0
            tx_npred = 0.0
            for line in fh:
                g = re.match(TEXON, line.strip())
                if g:
                    sen = 1.0 - float(g.group(1)) / int(g.group(2))
                    exon_sen.append(sen)

                g = re.match(FEXON, line.strip())
                if g:
                    tpp = 1.0 - float(g.group(1)) / int(g.group(2))
                    exon_tpp.append(tpp)

                g = re.match(TINTRON, line.strip())
                if g:
                    sen = 1.0 - float(g.group(1)) / int(g.group(2))
                    intron_sen.append(sen)

                g = re.match(FINTRON, line.strip())
                if g:
                    tpp = 1.0 - float(g.group(1)) / int(g.group(2))
                    intron_tpp.append(tpp)

                g = re.match(TLOCI, line.strip())
                if g:
                    sen = 1.0 - float(g.group(1)) / int(g.group(2))
                    loci_sen.append(sen)

                g = re.match(FLOCI, line.strip())
                if g:
                    tpp = 1.0 - float(g.group(1)) / int(g.group(2))
                    loci_tpp.append(tpp)

                g = re.match(TTX, line.strip())
                if g:
                    tx_tp = float(g.group(1))

                g = re.match(QTX, line.strip())
                if g:
                    tx_npred = float(g.group(1))

                g = re.match(RTX, line.strip())
                if g:
                    tx_nref = float(g.group(1))
            tx_sen.append(tx_tp / tx_nref)
            tx_tpp.append(tx_tp / tx_npred)

    print ("EXON_SEN\t{}".format("\t".join([str(i) for i in exon_sen])))
    print ("EXON_TPP\t{}".format("\t".join([str(i) for i in exon_tpp])))
    print ("INTRON_SEN\t{}".format("\t".join([str(i) for i in intron_sen])))
    print ("INTRON_TPP\t{}".format("\t".join([str(i) for i in intron_tpp])))
    print ("TX_SEN\t{}".format("\t".join([str(i) for i in tx_sen])))
    print ("TX_TPP\t{}".format("\t".join([str(i) for i in tx_tpp])))
    print ("LOCI_SEN\t{}".format("\t".join([str(i) for i in loci_sen])))
    print ("LOCI_TPP\t{}".format("\t".join([str(i) for i in loci_tpp])))
