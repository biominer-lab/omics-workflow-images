#!/usr/bin/env python
## created: 150202
## by: Joachim Weischenfeldt


# python /g/korbel//weischen/Dropbox/git/variant-calling/dellyVcf2Bedpe.py <vcf>

import os,sys,re,vcf, gzip, csv
from collections import defaultdict

if len(sys.argv) < 2:
	print 'syntax: <0> <vcf>'
	sys.exit(-1)

vcfFile = sys.argv[1]

if not 'vcf' in vcfFile:
	print 'file must contain vcf* suffix'
	sys.exit( -1)

#bedpeFile = vcfFile.split('.vcf')[0] +".bedpe"
bedpeFile = sys.argv[2]

print "input:", vcfFile

center = 'embl'

header = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "qual", "strand1", "strand2", "sv_class", "sv_type","score", "supp_reads","scna", "center", "read_id"]
with open(bedpeFile, 'w') as wout:
	writer = csv.writer(wout, delimiter="\t", lineterminator="\n")
	writer.writerow(header)
	reader = vcf.Reader(open(vcfFile, 'r'), compressed=True)
	sv_dict = defaultdict(dict)	
	for record in reader:
		if record.ID.endswith('_1'):
			coord1 = [record.CHROM, record.POS, int(record.POS)+1]
			if record.INFO.has_key('MATECHROM'):
				coord2 = [record.INFO['MATECHROM'], int(record.INFO['MATEPOS']), int(record.INFO['MATEPOS'])+1]
				strand2 = record.INFO['MATESTRAND']
			if record.INFO.has_key('CHR2'): # old format
				coord2 = [record.INFO['CHR2'], int(record.INFO['END']), int(record.INFO['END'])+1]
				strand2 = '+'
				if record.INFO['CT'].endswith("5"):
					strand2 = '-'
			svid = re.sub(r'_[12]$', '', record.ID)
			strand1 = record.INFO['STRAND']
			strands = [strand1, strand2]
			if 'SCNA_SUPP' in record.INFO.keys():
				scna_supp = '+'
			else:
				scna_supp = '-'
			supp_reads = record.INFO['PE']
			svinfo = [record.INFO['SVTYPE'], "BKND", "NA", supp_reads, scna_supp, center]
			try:
				readid = [','.join(record.INFO['READ_ID'])]
			except Exception, e:
				readid = ['NA']
			bedpe  = coord1 +  coord2 +  [svid, record.INFO['MAPQ']] + strands + svinfo  + readid
			writer.writerow(bedpe)

print "\n\nwrote to", bedpeFile

