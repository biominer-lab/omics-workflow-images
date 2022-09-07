#!/usr/bin/env python
## created: 150202
## by: Joachim Weischenfeldt

# python smufinVcf2Bedpe.py <vcf>

import os,sys,re,vcf, gzip, csv

from collections import defaultdict

if len(sys.argv) < 2:
    print 'syntax: <0> <vcf>'
    sys.exit(-1)

vcfFile = sys.argv[1]

if not 'vcf' in vcfFile:
    print 'file must contain vcf* suffix'
    sys.exit( -1)

# bedpeFile = vcfFile.split('.vcf')[0] +".bedpe"
bedpeFile = sys.argv[2]
print "input:", vcfFile

center = 'smufin'


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
        svid = re.sub(r'_[12]$', '', record.ID)
        strand1, strand2 = read_orient(record)
	#print (record.ALT, strand1, strand2)
        scna_supp = 'NA'
        qual = "NA"
        if record.QUAL:
            qual = record.QUAL
        if record.INFO.has_key('CNCH'):
            scna_supp = record.INFO['CNCH']
        readid = ['']
        if record.INFO.has_key('TSRDS'):
            readid = record.INFO['TSRDS']
        try:
            supp_reads = len(readid)
        except Exception, e:
            supp_reads = "NA"
        # svtype = record.INFO['SVTYPE']
        svinfo = ['NA', record.INFO['SVTYPE'], "NA", supp_reads, scna_supp, center] + [','.join(map(str, readid))]
        #svinfo = [record.INFO['SVTYPE'], ';'.join(map(str,record.INFO['SVCLASS'])), "NA", supp_reads, scna_supp, center]
        #svid = record.ID
        if record.ID.endswith('_1'):
            coord1 = [record.CHROM, record.POS, int(record.POS)+1]
            # coord1, qual, strand1, svinfo, readid = process_sanger(record)
            sv_dict[svid]["qual"] = qual
            sv_dict[svid]["coord1"] = coord1
            sv_dict[svid]["strand1"] = strand1
            sv_dict[svid]["svinfo"] = svinfo
            # sv_dict[svid]["readid"] = readid
        if record.ID.endswith('_2'):
            coord2 = [record.CHROM, record.POS, int(record.POS)+1]
            sv_dict[svid]["coord2"] = coord2
            sv_dict[svid]["strand2"] = strand1	
    for svkey in sv_dict.keys():
        bedpe = sv_dict[svkey]['coord1'] +  sv_dict[svkey]['coord2'] + [svkey, sv_dict[svkey]['qual'], sv_dict[svkey]['strand1'], sv_dict[svkey]['strand2']] + sv_dict[svkey]['svinfo'] 
        writer.writerow(bedpe)


print "\n\nwrote to", bedpeFile


