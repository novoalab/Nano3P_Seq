#!/usr/bin/env python
import sys
import pysam 
import os
import numpy as np 


'''
Sep. 25, 2021
Suzhou, China
yell @ huanle.liu@crg.eu if bugs are found
'''


usage = '''
	python clipped_read_ends_base_contents.py bamfile > outputfile
'''

if len (sys.argv) !=2:
	print (usage)
	sys.exit(1)

bam = sys.argv[1]

if not ( os.path.exists(bam) and os.path.isfile(bam)):
	print ('bam file does not exist; please ensure a correct input!')

if not os.path.isfile (bam+'.bai'):
	pysam.index (bam)

aln = pysam.AlignmentFile (bam, 'rb')

reflen = dict()
for x in aln.references:
	reflen[x] = aln.get_reference_length(x) 


print ('#ref\tread\tstrand\tsoftclip3end\tsoftclip5end\tsc3ATCG\tsc5ATCG')

for ref in aln.references:
	for rd in aln.fetch (contig = ref, start=0, stop=reflen[ref]):
		if rd.cigartuples is None or rd.query_sequence is None:
			continue  
		#if not rd.qname == '90a41d7a-052e-4ffc-a97e-ad4b273bcd9b':
		#	continue
		#else:
		#	print (rd.cigartuples)
		#	print (rd.get_aligned_pairs())
		#print (rd.get_aligned_pairs())	

		qseq = rd.query_sequence.upper()
		qlen = rd.query_length
		strand = '-' if rd.is_reverse else '+'
		length3, length5 = 0,0 
		op3, op5 = rd.cigartuples[0][0], rd.cigartuples[-1][0]
		
		if op3 == 4:
			length3 = rd.cigartuples[0][1]
		elif op5 == 4:
			length5 = rd.cigartuples[-1][1]

		clip3, clip5 = 'NONE','NONE'
		if strand == '-':
			length3, length5 = length5, length3
			op3, op5 = op5, op3
			if op3 == 4:
				clip3 = qseq[qlen - length3:]
			if op5 == 4:
				clip5 = qseq[:length5]
		elif strand == '+':
			if op3 == 4:
				clip3 = qseq[:length3]
			if op5 == 4:
				clip5 = qseq[qlen-length5:]
		clip3baseCNT = "-".join (map (str, [clip3.count(x) for x in 'ATCG'])) if clip3 != 'NONE' else 'NONE'
		clip5baseCNT = "-".join (map (str, [clip5.count(x) for x in 'ATCG'])) if clip5 != 'NONE' else 'NONE'
	
		print (rd.reference_name, rd.qname, strand, clip3, 
			clip5, clip3baseCNT, clip5baseCNT, sep='\t')
	






