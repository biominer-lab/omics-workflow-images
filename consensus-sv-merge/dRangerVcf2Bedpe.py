#!/usr/bin/env python
## created: 150202
## by: Joachim Weischenfeldt


# python dRangerVcf2Bedpe.py <vcf>
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
center = 'dRanger'

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
        if record.ID.endswith(':1'):
            coord1 = [record.CHROM, record.POS, int(record.POS)+1]
            matechrom = record.CHROM
            if ['MATECHROM'] in record.INFO.keys():
                matechrom = record.INFO['MATECHROM']
            elif record.ALT:
                #matechrom = re.sub(r'.*([^\:]):([0-9]*).*', r'\1' ,str(record.ALT))
                matechrom = re.sub(r'[^0-9XYZM]*([^:]*).*', r'\1', str(record.ALT))
            if ['MATEPOS'] in record.INFO.keys():
                pos2 = int(record.INFO['MATEPOS'])
            else:
                pos2 = int(re.sub(r'.*:([0-9]*).*', r'\1' ,str(record.ALT)))
            coord2 = [matechrom, int(pos2), int(pos2)+1]
            svid = re.sub(r':[12]$', '', record.ID)
            if ['STRAND'] in record.INFO.keys():
                strands = [record.INFO['STRAND'], record.INFO['MATESTRAND']]
            else:
                strands = read_orient(record)    
            scna_supp = 'NA'
            supp_reads = 'NA'
            svclass = 'NA'
            if ['SVCLASS'] in record.INFO.keys():
                svclass = record.INFO['SVCLASS']
            somaticScore = "NA"
            if ['SOMATICSCORE'] in record.INFO.keys():
                somaticScore = record.INFO['SOMATICSCORE']
            qual = "NA"
            if ['DRQUAL'] in record.INFO.keys():
                qual = record.INFO['DRQUAL']
            elif ['MAPQ'] in record.INFO.keys():
                qual = record.INFO['MAPQ']
            readid = ['.']
            try:
                readid = [','.join([call['READ_ID'] for call in record.samples if call['READ_ID']][0])]
                supp_reads = len(readid[0].split(','))
            except Exception, e:
                pass
            svinfo = [svclass, record.INFO['SVTYPE'], somaticScore, supp_reads, scna_supp, center]
            bedpe  = coord1 +  coord2 +  [svid, qual] + strands + svinfo  + readid 
            writer.writerow(bedpe)

print ("\n\nwrote to", bedpeFile)
