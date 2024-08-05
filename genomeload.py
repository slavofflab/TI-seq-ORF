"""
Script to load genome fa file.
"""

from Bio import SeqIO
from Bio.Seq import Seq

def Genomeload(gene_file):
    if  gene_file.split(".")[-1] == "gz":
        try:
            gf = gzip.open(gene_file,'r')
        except IOError:
            print("The gene file does not exist!")
            parser.print_help()
            sys.exit()
    else:
        try:
            gf = open(gene_file,'r')
        except IOError:
            print("The gene file does not exist!")
            parser.print_help()
            sys.exit()

    chrom_seq = {}

    print("Start loading gene file")
    for record in SeqIO.parse(gf,"fasta"):
        chrom_seq[record.id] = Seq(str(record.seq)).upper()
        print(record.id + " loaded")

    gf.close()
    return chrom_seq
