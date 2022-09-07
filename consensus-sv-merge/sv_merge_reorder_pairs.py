#!/usr/bin/env python
## created: 150204
## by: Joachim Weischenfeldt

#python sv_merging/pcawg_merge_reorder_pairs.py  <file> <pid>


from __future__ import division
import os,sys,re,vcf, gzip, csv
from collections import defaultdict

merge_file = sys.argv[1]
pid = sys.argv[2]

#print "reorder merged SV calls"
slop=500.0
header = ["chrom1_SV1", "pos1_SV1", "chrom2_SV1", "pos2_SV1", "name_SV1", "score_SV1", "strand1_SV1", "strand2_SV1", "svtype_SV1", "center_SV1", "idSV1", "chrom1_SV2", "pos1_SV2", "chrom2_SV2", "pos2_SV2", "name_SV2", "score_SV2", "strand1_SV2", "strand2_SV2", "svtype_SV2", "center_SV2", "idSV2", "pid", "interSect", "BpOffset", "svOverlapType", "callerPair", "rn_share_count", "rn_share_fract"]

print ('\t'.join(map(str, header)))

## line size limit might be compromised. Increase
old_limit = csv.field_size_limit()
csv.field_size_limit(13107200)



with open(merge_file, "r") as rin:
    reader = csv.reader(rin, delimiter="\t", quoting=csv.QUOTE_NONE)
    center_seen = set()
    rowadj = 0
    for row in reader:
        if len(row) == 33:
            rowadj = -1
        if row[-1] == pid:
            center1 = row[15]
            center2 = row[32 + rowadj]
            centerset = '__'.join(map(str, sorted([row[6],row[23 + rowadj]])))
            if center1 != center2 and centerset not in center_seen: # not same caller
                center_seen.add(centerset)
                if row[0] == row[17 + rowadj]:
                    distPos1 = abs(int(row[1]) - int(row[18 + rowadj]))
                elif row[0] == row[20 + rowadj]:
                    distPos1 = abs(int(row[1]) - int(row[21 + rowadj]))
                if row[3] == row[20 + rowadj]:
                    distPos2 = abs(int(row[4]) - int(row[21 + rowadj]))
                elif row[3] == row[17]:
                    distPos2 = abs(int(row[4]) - int(row[18 + rowadj]))
                interSect = (4*slop - distPos1 - distPos2)/(4*slop)
                BpOffset = distPos1 + distPos2
                svType1 = row[10]
                if svType1 == "NA":
                    svType1 = row[11]
                svType2 = row[27]
                if svType2 == "NA":
                    svType2 = row[28 + rowadj]
                svOverlapType = svType1 + "_" + svType2
                callerPair = center1 + "_" +  svType1 + "|" + center2 + "_" + svType2
                idSV1 = center1 + "_" + row[6] + "_" + row[0] + ":" + row[1] + "-" + row[3] + ":" + row[4]
                idSV2 = center2 + "_" + row[23+ rowadj] + "_" + row[17] + ":" + row[18 + rowadj] + "-" + row[20+ rowadj] + ":" + row[21+ rowadj]
                rn_SV1 = set([i.strip('/[12]').replace('-', ':') for i in row[16].split(',')])
                rn_SV2 = set([i.strip('/[12]').replace('-', ':') for i in row[33 + rowadj].split(',')])
                rn_union = rn_SV1.union(rn_SV2)
                rn_shared = rn_SV1.intersection(rn_SV2)
                rn_share_count = len(rn_shared)
                if not 'NA' in rn_SV1 and not 'NA' in rn_SV2 and not rn_shared:
                    continue
                    # kick it out if each caller has readid but none are shared
                try:
                    rn_share_fract = round(rn_share_count/len(rn_union),2)
                except Exception, e:
                    rn_share_fract = 0
                if rowadj == -1:
                    out = row[0:2] + row[3:5] + row[6:11] + [center1, idSV1] + row[16:18] + row[19:21] + row[22:27] + [center2, idSV2, pid] + [interSect, BpOffset, svOverlapType, callerPair, rn_share_count, rn_share_fract]
                else:
                    out = row[0:2] + row[3:5] + row[6:11] + [center1, idSV1] + row[17:19] + row[20:22] + row[23:28] + [center2, idSV2, pid] + [interSect, BpOffset, svOverlapType, callerPair, rn_share_count, rn_share_fract]
                print ('\t'.join(map(str, out)))
