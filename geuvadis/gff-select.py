#!/home/ruolin/Software/anaconda3/bin/python
import sys
import argparse
import re

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("gff", type=str, help = "the gff file")
    parser.add_argument("-l","--list", type=str, required=True, help="a list of gff name to be selected")
    parser.add_argument("-g","--gene", default=False, action='store_true', help="provide a list of gene id. Default is a transcript id")
    args = parser.parse_args()

    gff_f = args.gff
    subset_f = args.list

    gene = re.compile("gene_id \"(.*?)\";")
    iso = re.compile("transcript_id=(.*?);")

    chosen = {}

    with open(subset_f, 'r') as fh:
        for line in fh:
            chosen[line.strip()] = 1

    with open(gff_f, 'r') as fh:
        for line in fh:
            if line.startswith('#'): continue
            else:
                fields = line.strip().split("\t")
                ret = None
                if args.gene:
                    ret = gene.search(fields[8])
                else:
                    ret = iso.search(fields[8])
                if ret:
                    t = ret.group(1) 
                    if t in chosen:
                        print(line.strip())
