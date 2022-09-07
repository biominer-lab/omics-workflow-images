#!/usr/bin/env python
## created: 220902
## by: Jun Shang


# python mantaVcf2Bedpe.py <vcf>
from __future__ import print_function
import os,sys,re,vcf, gzip, csv
from collections import defaultdict

if len(sys.argv) < 2:
    print ('syntax: <0> <vcf>')
    sys.exit(-1)

vcfFile = sys.argv[1]

if not 'vcf' in vcfFile:
    print ('file must contain vcf* suffix')
    sys.exit( -1)

#bedpeFile = vcfFile.split('.vcf')[0] +".sv.bedpe"
bedpeFile = sys.argv[2] 

print ("input:", vcfFile)

center = vcfFile.split('/')[-2]
center = 'manta'

def read_orient(record):
    strands = ['NA', 'NA']
    if re.search(r'^[ACTGN]+\]',str(record.ALT[0])):
        strands = ['+', '+']
    if re.search(r'^[ACTGN]+\[',str(record.ALT[0])):
        strands = ['+', '-']
    if re.search(r'\][ACTGN]+$',str(record.ALT[0])):
        strands = ['-', '+']
    if re.search(r'\[[ACTGN]+$',str(record.ALT[0])):
        strands = ['-', '-']
    return strands


header = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "qual", "strand1", "strand2", "sv_class", "sv_type","score", "supp_reads","scna", "center", "read_id"]
with open(bedpeFile, 'w') as wout:
    writer = csv.writer(wout, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    reader = vcf.Reader(open(vcfFile), 'r', compressed=True) if vcfFile.endswith('.gz') else vcf.Reader(open(vcfFile), 'r', compressed=False)
    sv_dict = defaultdict(dict)    
    for record in reader:
        svid = re.sub(r':[01]$', '', record.ID)
        strand1, strand2 = read_orient(record)
        qual = "NA"
        sv_class = "NA"
        score = "NA"
        scna_supp = "NA"
        supp_reads = "NA"
        read_id = "NA"
        
        if ['SVTYPE'] in record.INFO.keys():
            sv_type = record.INFO['SVTYPE']
        
        svinfo = [sv_class, record.INFO['SVTYPE'], score, supp_reads, scna_supp, center, read_id]
        
        if record.ID.endswith(':0'):
            coord1 = [record.CHROM, record.POS, int(record.POS)+1]
            # coord1, qual, strand1, svinfo, readid
            sv_dict[svid]["qual"] = qual
            sv_dict[svid]["coord1"] = coord1
            sv_dict[svid]["strand1"] = strand1
            sv_dict[svid]["svinfo"] = svinfo
            
        if record.ID.endswith(':1'):
            coord2 = [record.CHROM, record.POS, int(record.POS)+1]
            sv_dict[svid]["coord2"] = coord2
            sv_dict[svid]["strand2"] = strand1
            
    for svkey in sv_dict.keys():
        if all(k in sv_dict[svkey] for k in ("coord1","coord2")):
            bedpe = sv_dict[svkey]['coord1'] +  sv_dict[svkey]['coord2'] + [svkey, sv_dict[svkey]['qual'], sv_dict[svkey]['strand1'], sv_dict[svkey]['strand2']] + sv_dict[svkey]['svinfo']
        else:
            continue
        writer.writerow(bedpe)
    
print("\n\nwrote to", bedpeFile)
