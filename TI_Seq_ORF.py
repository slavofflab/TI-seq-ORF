from genomeload import Genomeload
from gtfline import gtf_line
from mergeinterval import merge
import re
import warnings

warnings.filterwarnings('ignore')

def maxlocal(RL):
    M = RL[3]
    for k in range(0,7):
        if k!= 3 and RL[k] >= M:
            return False
    return True

print("Start loading Genome")
genome = Genomeload("pri_hg38.fa")
readlist = {}
for chrom in genome:
    readlist[chrom + "+"] = [0.0]*len(genome[chrom])
    readlist[chrom + "-"] = readlist[chrom + "+"]
print("Finish loading Genome")

print("Start loading TI-seq reads")
with open("HEK293_LTM1.+.bedgraph", "r") as readfile1:
    for line in readfile1:
        temp = line.strip().split("\t")
        for i in range(int(temp[1])+1,int(temp[1])+2):
            readlist["chr" + temp[0] + "+"][i] = float(temp[3])

with open("HEK293_LTM1.-.bedgraph", "r") as readfile2:
    for line in readfile2:
        temp = line.strip().split("\t")
        for i in range(int(temp[1])+3,int(temp[1])+4):
            readlist["chr" + temp[0] + "-"][i] = float(temp[3])
print("Finish loading TI-seq reads")

print("Start loading gtf file")
rna_set = {}
annotation_file = "gencode.v38.annotation.gtf"

for key in genome:
    rna_set[key + "+"] = {}
    rna_set[key + "-"] = {}

if annotation_file.split(".")[-1] == "gz":
    try:
        af = gzip.open(annotation_file,'r')
    except IOError:
        print("The annotation file does not exist!")
        parser.print_help()
        sys.exit()
else:
    try:
        af = open(annotation_file,'r')
    except IOError:
        print("The annotation file does not exist!")
        parser.print_help()
        sys.exit()

content = af.readline()
tc = 0

while content:
    if len(content.split("\t")) < 8:
        content = af.readline()
    else:
        gtf_list = gtf_line(content.strip())
        if gtf_list.type == "gene":
            rna_set[gtf_list.chrom + gtf_list.strand][gtf_list.attributes["gene_name"]] = {}
            content = af.readline()
        elif gtf_list.type == "transcript":
            rna_set[gtf_list.chrom + gtf_list.strand][gtf_list.attributes["gene_name"]][gtf_list.attributes["transcript_id"]] = {"TSS":gtf_list.start, "TES":gtf_list.end, "exons":[], "CDS":[]}
            content = af.readline()
            tc += 1
            if tc%10000 == 0:
                print(f"Loaded {tc} transcripts")
        elif gtf_list.type == "exon":
            rna_set[gtf_list.chrom + gtf_list.strand][gtf_list.attributes["gene_name"]][gtf_list.attributes["transcript_id"]]["exons"].append([gtf_list.start, gtf_list.end])
            content = af.readline()
        elif gtf_list.type == "CDS":
            temCDS = rna_set[gtf_list.chrom + gtf_list.strand][gtf_list.attributes["gene_name"]][gtf_list.attributes["transcript_id"]]["CDS"]
            rna_set[gtf_list.chrom + gtf_list.strand][gtf_list.attributes["gene_name"]][gtf_list.attributes["transcript_id"]]["CDS"] = [min(temCDS + [gtf_list.start]), max(temCDS + [gtf_list.end])]
            content = af.readline()
        else:
            content = af.readline()
af.close()

print("Loaded " + str(tc) + " transcripts")
print("Start analyzing reads")
out1, out2 = open("ateORF.out", "w"), open("TIS.out", "w")
compatible = ["ncORF","uORF","uoORF","novel"]
out1.write(f"Chromosome\tStrand\tPosition\tIntensity\tGene name\tObstructed by CDS in all transcripts\tRange\tStart codon\tStart codon context\tKozak\tORF max length\tAll ORF types\tpeptides sequences\n")
out2.write(f"Chromosome\tStrand\tPosition\tIntensity\tGene name\tORF type\tRange\tStart codon\tStart codon context\tKozak\tORF max length\tAll ORF types\tpeptides sequences\n")

top,threshold = 5, 5
for ch, ch_set in rna_set.items():
    for gene, gene_set in ch_set.items():
        exon_merged, start = [], []
        for trans, trans_set in gene_set.items():
            #Merge exon regions
            exon_merged += [[v[0], v[1]] for v in trans_set["exons"]]
            exon_merged = merge(exon_merged)
        posLIST,readLIST = [], []
        for exon in exon_merged:
            posLIST += range(exon[0], exon[1] + 1)
            readLIST += readlist[ch][exon[0]:exon[1] + 1]
        ziplist = list(zip(readLIST, posLIST))
        ziplist.sort(reverse=True)
        i = 0
        h = ziplist[0][0]
        while i <= top:
            if ziplist[i][0] <= threshold or ziplist[i][0] <= 0.2*h:
                #Filting low intensity peaks
                break
            else:
                temlist = ([0,0,0]+readLIST+[0,0,0])[posLIST.index(ziplist[i][1]):posLIST.index(ziplist[i][1])+7]
                if maxlocal(temlist) == False:
                    i += 1
                    pass
                else:
                    peak = ziplist[i][1]
                    typelist, peplist, ORFmaxlength = [],[], 0
                    sc0, sc1 = "", ""
                    for trans, trans_set in gene_set.items():
                        CDS = trans_set["CDS"]
                        exons = trans_set["exons"]
                        exons.sort()
                        POSLIST, seqLIST = [], ""
                        for exon in exons:
                            POSLIST += range(exon[0], exon[1] + 1)
                            seqLIST += genome[ch[:-1]][exon[0]-1:exon[1]]
                        if ch[-1] == "-":
                            POSLIST = POSLIST[::-1]
                            seqLIST = seqLIST.reverse_complement()
                            CDS = CDS[::-1]
                        if peak not in POSLIST:
                            pass
                        else:
                            POSslice = POSLIST[POSLIST.index(peak):]
                            SEQslice = seqLIST[POSLIST.index(peak):]
                            sc00 = SEQslice[0:3]
                            if len(sc0) < len(sc00):
                                sc0 = sc00
                            sc10 = seqLIST[POSLIST.index(peak)-3:POSLIST.index(peak)+4]
                            if len(sc1) < len(sc10):
                                sc1 = sc10
                            Kozak = "No context"
                            if re.match(r'[AG]..ATGG', str(sc1)):
                                Kozak = "Strong"
                            elif re.match(r'[CT]..ATGG', str(sc1)) or re.match(r'[AG]..ATG[ACT]', str(sc1)):
                                Kozak = "Adequate"
                            elif re.match(r'[CT]..ATG[ACT]', str(sc1)):
                                Kozak = "Weak"
                            pep = SEQslice.translate()
                            if "*" in pep:
                                pep = pep.split("*")[0]
                                if pep not in peplist:
                                    peplist.append(str(pep))
                                    if len(pep) > ORFmaxlength:
                                        ORFmaxlength = len(pep)
                                if CDS == []:
                                    typelist.append("ncORF")
                                elif POSLIST.index(peak) >  POSLIST.index(CDS[1]):
                                    typelist.append("dORF")
                                elif POSLIST.index(peak) == POSLIST.index(CDS[0]):
                                    typelist.append("CDS")
                                elif POSLIST.index(peak) < POSLIST.index(CDS[0]):
                                    sp = POSLIST.index(POSslice[len(pep)*3])
                                    if sp == (POSLIST.index(CDS[1])+1):
                                        typelist.append("N-ext")
                                    elif sp < POSLIST.index(CDS[0]):
                                        typelist.append("uORF")
                                    elif sp < (POSLIST.index(CDS[1])+1):
                                        typelist.append("uoORF")
                                    else:
                                        typelist.append("novel")
                                elif POSLIST.index(peak) < POSLIST.index(CDS[1]):
                                    sp = POSLIST.index(POSslice[len(pep)*3])
                                    if sp == POSLIST.index(CDS[1])+1:
                                        typelist.append("N-del")
                                    else:
                                        typelist.append("iORF")
                            else:
                                pass
                                #No stop ORF
                    if "CDS" in typelist:
                        ORFtype = "CDS"
                    elif "N-ext" in typelist or "N-del" in typelist:
                        ORFtype = "isoform"
                    elif "iORF" in typelist:
                        ORFtype = "iORF"
                        if ORFmaxlength > 5:
                            if any(ot in typelist for ot in compatible):
                                Obstruction  = "NO"
                            else:
                                Obstruction  = "Yes"
                            out1.write(f"{ch[:-1]}\t{ch[-1]}\t{ziplist[i][1]}\t{ziplist[i][0]}\t{gene}\t{Obstruction}\t{ch[:-1]}:{ziplist[i][1]-5}-{ziplist[i][1]+5}\t{sc0}\t{sc1}\t{Kozak}\t{ORFmaxlength}\t{typelist}\t{peplist}\n")
                    elif "dORF" in typelist:
                        ORFtype = "dORF"
                    elif "uORF" in typelist:
                        ORFtype = "uORF"
                    elif "uoORF" in typelist:
                        ORFtype = "uoORF"
                    elif "novel" in typelist:
                        ORFtype = "novel"
                    elif "ncORF" in typelist:
                        ORFtype = "ncORF"
                    else:
                        ORFtype = "no transcript"
                    if ORFmaxlength > 5:
                        out2.write(f"{ch[:-1]}\t{ch[-1]}\t{ziplist[i][1]}\t{ziplist[i][0]}\t{gene}\t{ORFtype}\t{ch[:-1]}:{ziplist[i][1]-5}-{ziplist[i][1]+5}\t{sc0}\t{sc1}\t{Kozak}\t{ORFmaxlength}\t{typelist}\t{peplist}\n")
                    i += 1 
out1.close()
out2.close()
print("Finished analyzing reads")
