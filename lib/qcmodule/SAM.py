#!/usr/bin/env python
'''manipulate BAM/SAM file.'''

#import built-in modules
import os,sys
import re
import string
from optparse import OptionParser
import warnings
import string
import collections
import math
import sets
import random

#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *
from bx.binned_array import BinnedArray
from bx_extras.fpconst import isNaN
from bx.bitset_utils import *
import pysam
from qcmodule import mystat
from qcmodule import fasta
from qcmodule import bam_cigar
from qcmodule import BED
from qcmodule import Results

#changes to the paths

class ReadDup_prof:
	def __init__(self,outfile_data,outfile_fig):
		self.outfile1 = outfile_data+".seq.DupRate.xls"
		self.outfile2 = outfile_data + ".pos.DupRate.xls"
		self.outfile3 = outfile_data + ".DupRate_plot.r"
		
		self.seqDup=collections.defaultdict(int)
		self.posDup=collections.defaultdict(int)
		
		self.plot_file = outfile_fig+".DupRate_plot.png"
		
		
	def write(self,up_bound=500):
		'''Calculate reads's duplicate rates'''
		SEQ=open(self.outfile1,'w')
		POS=open(self.outfile2,'w')
		RS=open(self.outfile3,'w')

		seqDup_count=collections.defaultdict(int)
		posDup_count=collections.defaultdict(int)		

		print >>SEQ, "Occurrence\tUniqReadNumber"
		for i in self.seqDup.values():			#key is occurence, value is uniq reads number (based on seq)
			seqDup_count[i] +=1
		for k in sorted(seqDup_count.iterkeys()):	
			print >>SEQ, str(k) +'\t'+ str(seqDup_count[k])
		SEQ.close()
		
		print >>POS, "Occurrence\tUniqReadNumber"
		for i in self.posDup.values():			#key is occurence, value is uniq reads number (based on coord)
			posDup_count[i] +=1
		for k in sorted(posDup_count.iterkeys()):	
			print >>POS, str(k) +'\t'+ str(posDup_count[k])
		POS.close()
		
		print >>RS, "png(\'%s\',width=600,height=600,units='px')" % (self.plot_file)
		print >>RS, "par(mar=c(5,4,4,5),las=0)"
		print >>RS, "seq_occ=c(" + ','.join([str(i) for i in sorted(seqDup_count.iterkeys()) ]) + ')'
		print >>RS, "seq_uniqRead=c(" + ','.join([str(seqDup_count[i]) for i in sorted(seqDup_count.iterkeys()) ]) + ')'
		print >>RS, "pos_occ=c(" + ','.join([str(i) for i in sorted(posDup_count.iterkeys()) ]) + ')'
		print >>RS, "pos_uniqRead=c(" + ','.join([str(posDup_count[i]) for i in sorted(posDup_count.iterkeys()) ]) + ')'
		print >>RS, "plot(pos_occ,log10(pos_uniqRead),ylab='Number of Reads (log10)',xlab='Frequency',pch=4,cex=0.8,col='blue',xlim=c(1,%d),yaxt='n')" % up_bound
		print >>RS, "points(seq_occ,log10(seq_uniqRead),pch=20,cex=0.8,col='red')"
		print >>RS, 'ym=floor(max(log10(pos_uniqRead)))'
		print >>RS, "legend(%d,ym,legend=c('Sequence-base','Mapping-base'),col=c('blue','red'),pch=c(4,20))" % max(up_bound-200,1)
		print >>RS, 'axis(side=2,at=0:ym,labels=0:ym)'
		print >>RS, 'axis(side=4,at=c(log10(pos_uniqRead[1]),log10(pos_uniqRead[2]),log10(pos_uniqRead[3]),log10(pos_uniqRead[4])), labels=c(round(pos_uniqRead[1]*100/sum(pos_uniqRead)),round(pos_uniqRead[2]*100/sum(pos_uniqRead)),round(pos_uniqRead[3]*100/sum(pos_uniqRead)),round(pos_uniqRead[4]*100/sum(pos_uniqRead))))'
		print >>RS, 'mtext(4, text = "Reads %", line = 2)'
		print >>RS, 'dev.state = dev.off()'
		
		RS.close()
	
	def count(self,RNA_read,chrom,hit_st,exon_blocks):
				exon_boundary=""
						
				self.seqDup[RNA_read] +=1					#key is read sequence				
				for ex in exon_blocks:
					exon_boundary += str(ex[1]) + '-' + str(ex[2]) + ":"
				key = chrom + ":" + str(hit_st) + ":" + exon_boundary
				self.posDup[key] +=1


class InnerDist_prof:
	def __init__(self,outfile_data,outfile_fig,sample_size=1000000,low_bound=0,up_bound=1000,step=10):
		
		self.out_file1 = outfile_data + ".inner_distance.txt"	
		self.out_file2 = outfile_data + ".inner_distance_freq.txt"
		self.out_file3 = outfile_data + ".inner_distance_plot.r"
		self.plot_file = outfile_fig+".inner_distance_plot.png"
		
		self.ranges ={}
		self.ranges['chr100'] = Intersecter()
		self.pair_num = 0
		self.sample_size = sample_size
		self.step = step
		self.window_left_bound = range(low_bound,up_bound,step)		
		
	def write(self):
		'''estimate the inner distance of mRNA pair end fragment. fragment size = insert_size + 2 x read_length'''
		if self.pair_num == 0 :
			return
			
		try :
			#FO=open(self.out_file1,'w')
			FQ=open(self.out_file2,'w')
			RS=open(self.out_file3,'w')
		except:
			sys.stderr.write("Error in creating file %s \n" %self.out_file1)
			pass
		
		sizes = []
		counts = []
		
		for st in self.window_left_bound:
			sizes.append(str(st + self.step/2))
			count = str(len(self.ranges['chr100'].find(st,st + self.step)))
			counts.append(count)
			print >>FQ, str(st) + '\t' + str(st+self.step) +'\t' + count		
		
		#plot_file = outfile_fig+".inner_distance_plot.png"
		print >>RS, "png(\"%s\",width=600,height=600,units='px')" % (self.plot_file)
		#print >>RS, "par(mfrow=c(2,1),cex.main=0.8,cex.lab=0.8,cex.axis=0.8,mar=c(4,4,4,1))"
		#print >>RS, 'pie(c(%d,%d,%d),col=rainbow(3),cex=0.5,radius=1,main="Total %d fragments",labels=c("fraSize <= %d\\n(%4.2f%%)","fragSize > %d\\n(%4.2f%%)","%d < fragSize <= %d\\n(%4.2f%%)"), density=rep(80,80,80),angle=c(90,140,170))' % (ultra_low, ultra_high, pair_num -ultra_low -ultra_high, pair_num, low_bound, ultra_low*100/pair_num, up_bound, ultra_high*100/pair_num, low_bound, up_bound, 100-ultra_low*100/pair_num - ultra_high*100/pair_num)
		print >>RS, 'fragsize=rep(c(' + ','.join(sizes) + '),' + 'times=c(' + ','.join(counts) + '))'
		print >>RS, 'frag_sd = sd(fragsize)'
		print >>RS, 'frag_mean = mean(fragsize)'
		print >>RS, 'frag_median = median(fragsize)'
		print >>RS, 'write(c("Mean insert size",frag_mean), stdout())'
		print >>RS, 'write(c("Median insert size",frag_median), stdout())'
		print >>RS, 'write(c("Standard deviation",frag_sd), stdout())'
		print >>RS, 'hist(fragsize,probability=T,breaks=%d,xlab="mRNA insert size (bp)",main=paste(c("Mean=",frag_mean,";","SD=",frag_sd),collapse=""),border="blue")' % len(self.window_left_bound)
		print >>RS, "lines(density(fragsize,bw=%d),col='red')" % (2*self.step)
		print >>RS ,"dev.state = dev.off()"
		#FO.close()
		FQ.close()
		RS.close()
		
		
	def count(self,exon_bitsets,aligned_read,chrom):
		
		if self.pair_num > self.sample_size:
			return
		self.pair_num +=1
		
		frag_size=0
		fchrom = 'chr100'	
		
		inner_distance_bitsets=BinnedBitSet()
		tmp = BinnedBitSet()
		tmp.set_range(0,0)

		splice_intron_size=0			

		read1_len = aligned_read.qlen
		read1_start = aligned_read.pos
		read2_start = aligned_read.mpos
				
				
		#chrom = self.samfile.getrname(aligned_read.tid).upper()
		intron_blocks = bam_cigar.fetch_intron(chrom, read1_start, aligned_read.cigar)				
		for intron in intron_blocks:
				splice_intron_size += intron[2] - intron[1]
				
		read1_end = read1_start + read1_len + splice_intron_size		
				#if read1_end > read2_start: continue
				
		inner_distance = read2_start - read1_end +1
		
		if inner_distance > 0: 
					if chrom in exon_bitsets:
						size =0 
						inner_distance_bitsets.set_range(read1_end, read2_start-read1_end)
						inner_distance_bitsets.iand(exon_bitsets[chrom])
						end=0
						while 1:
							start = inner_distance_bitsets.next_set( end )
							if start == inner_distance_bitsets.size: break
							end = inner_distance_bitsets.next_clear( start )
							size += (end - start) +1
						inner_distance_bitsets.iand(tmp)											#clear BinnedBitSet
						if size == inner_distance:
							#FO.write(aligned_read.qname + '\t' + str(size) + '\tPE_within_same_exon\n')
							self.ranges[fchrom].add_interval( Interval( size-1, size ) )
						elif size < inner_distance and size >1:
							#if same_chrom:
								#FO.write(aligned_read.qname + '\t' + str(size) + '\tPE_within_diff_exon\n')
								self.ranges[fchrom].add_interval( Interval( size-1, size ) )	

						elif size == 1:
							#FO.write(aligned_read.qname + '\t' + str(inner_distance) + '\tPE_not_in_gene_region\n')
							self.ranges[fchrom].add_interval( Interval( inner_distance-1, inner_distance ) )
						
					else:
						#FO.write(aligned_read.qname + '\t' + str(inner_distance) + '\tPE_chrom_has_no_refgene\n')
						self.ranges[fchrom].add_interval( Interval( inner_distance-1, inner_distance ) )
		else:
					#FO.write(aligned_read.qname + '\t' + str(inner_distance) + '\tPE_reads_overlap\n')
					self.ranges[fchrom].add_interval( Interval( inner_distance-1, inner_distance ) )
	

class CoverageBody_prof:
	def __init__(self,outfile_data,outfile_fig,exon_r):
		self.outfile1 = outfile_data + ".geneBodyCoverage_plot.r"
		self.outfile2 = outfile_data + ".geneBodyCoverage.txt"
		self.plot_file = outfile_fig + ".geneBodyCoverage.png"
		self.outfile3 = outfile_data + ".geneAbundance.txt"
		
		self.ranges = {}
		self.frag_num = 0
		self.geneCounts_list = dict()
		self.annotations = exon_r
		
	def write(self,totalReads,refbed):

		coverage=collections.defaultdict(int)	
		#geneCounts_list = dict()
		
		for line in open(refbed,'r'):
			try:
				if line.startswith(('#','track','browser')):continue  
				# Parse fields from gene tabls
				fields = line.split()
				chrom	 = fields[0].upper()
				tx_start  = int( fields[1] )
				tx_end	= int( fields[2] )
				geneName	  = fields[3]
				strand	= fields[5]
				
				exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
				exon_starts = map((lambda x: x + tx_start ), exon_starts)
				exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
				exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends);   
			except:
				print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line,
				continue
			gene_all_base=[]
			percentile_base=[]
			mRNA_len =0			
			
			for st,end in zip(exon_starts,exon_ends):
				gene_all_base.extend(range(st+1,end+1))		#0-based coordinates on genome
				mRNA_len = len(gene_all_base)
			if mRNA_len <100:
				continue
			if strand == '-':
				gene_all_base.sort(reverse=True)			#deal with gene on minus stand
			else:
				gene_all_base.sort(reverse=False)
			percentile_base = mystat.percentile_list (gene_all_base)	#get 101 points from each gene's coordinates
			
			#hit = 0
			for i in range(0,len(percentile_base)):
				if chrom in self.ranges:
					cc =  len(self.ranges[chrom].find(percentile_base[i], percentile_base[i]+1))
					coverage[i] += cc	
					
		x_coord=[]
		y_coord=[]
		try :
			OUT2 = open(self.outfile2,'w')
			OUT1 = open(self.outfile1,'w')
			OUT3 = open(self.outfile3,'w')
			
			print >>OUT2, "Total reads: " + str(totalReads)
			print >>OUT2, "Fragment number: " + str(self.frag_num)
			print >>OUT2, "percentile\tcount"
			
			for i in range(len(coverage)):
				x_coord.append(str(i))
				y_coord.append(str(coverage[i]))
				print >>OUT2, str(i) + '\t' + str(coverage[i])
			print >>OUT1, "png(\'%s\',width=600,height=600,units='px')" % (self.plot_file)
			print >>OUT1, "x=c("+','.join(x_coord)+')'
			print >>OUT1, "y=c(" + ','.join(y_coord) + ')'
			print >>OUT1, "plot(x,y,xlab=\"percentile of gene body (5'->3')\",ylab='read number',type='s',xlim=c(0,100))"
			print >>OUT1, "dev.state = dev.off()"
			
			OUT1.close()
			OUT2.close()
			
			if len(self.geneCounts_list) >0 :
				#sys.stderr.write(str(len(self.geneCounts_list))+"\n")
				for key in self.geneCounts_list :
					#sys.stderr.write(key+"\n")
					print >>OUT3, key +"\t"+str(self.geneCounts_list[key])
				
			OUT3.close()
						
		except :
			sys.stderr.write("Error in writing gene body coverage data.\n")
			pass		

	def count2(self,exon_blocks,exon_blocks2,chrom1,chrom2):	
		#sys.stderr.write(str(w)+"\n")			
		self.frag_num += len(exon_blocks)+len(exon_blocks)
		genes1 = []
		genes2 = []
		for exon in exon_blocks:
			if chrom1 not in self.ranges:
				self.ranges[chrom1] = Intersecter()
			self.ranges[chrom1].add_interval( Interval( exon[1], exon[2] ) )
			if chrom1 in self.annotations:
				found_list = self.annotations[chrom1].find(exon[1],exon[2])
				for i in range(len(found_list)) :					
					geneName = found_list[i]
					#sys.stderr.write(geneName+"\n")
					if geneName not in genes1 :
						genes1.append(geneName)
		
		for exon in exon_blocks2:
			if chrom2 not in self.ranges:
				self.ranges[chrom2] = Intersecter()
			self.ranges[chrom2].add_interval( Interval( exon[1], exon[2] ) )
			if chrom2 in self.annotations:
				found_list = self.annotations[chrom2].find(exon[1],exon[2])
				for i in range(len(found_list)) :					
					geneName = found_list[i]
					#sys.stderr.write(geneName+"\n")
					if geneName not in genes2 :
						genes2.append(geneName)
						
		#intersect genes1 and genes2
		genes = []
		
		if len(genes1) ==0 or len(genes2) == 0 :
			genes.extend(genes1)
			genes.extend(genes2)
		else :
			genes = list(set(genes1).intersection(set(genes2)))	
		#for g in genes :
		if len(genes) == 1 : #no ambiguous assignment
			g = genes[0]
			if g not in self.geneCounts_list :
				self.geneCounts_list[g] = 1.0
			else :
				self.geneCounts_list[g] += 1.0
					
	def count(self,exon_blocks,chrom):	
		#sys.stderr.write(str(w)+"\n")			
		self.frag_num += len(exon_blocks)
		genes = []
		for exon in exon_blocks:
			if chrom not in self.ranges:
				self.ranges[chrom] = Intersecter()
			self.ranges[chrom].add_interval( Interval( exon[1], exon[2] ) )
			if chrom in self.annotations:
				found_list = self.annotations[chrom].find(exon[1],exon[2])
				for i in range(len(found_list)) :					
					geneName = found_list[i]
					#geneName = iv.value
					#sys.stderr.write(geneName+"\n")
					if geneName not in genes :
						genes.append(geneName)
		
		#for g in genes :
		if len(genes) == 1 : #no ambiguous assignment
			g = genes[0]
			if g not in self.geneCounts_list :
				self.geneCounts_list[g] = 1.0/len(genes)
			else :
				self.geneCounts_list[g] += 1.0/len(genes)
		
class Clipping_prof:
	def __init__(self,outfile_data,outfile_fig):
		
		self._out_file1 = outfile_data + ".clipping_profile.xls"
		self._out_file2 = outfile_data + ".clipping_profile.r"
		self.plot_file = outfile_fig+ ".clipping_profile.png"
		
		self.soft_p = re.compile(r'(.*?)(\d+)S')
		self.read_part = re.compile(r'(\d+)[MIS=X]')
		#total_read =0
		#self.skip_part_of_read =0
		self.soft_clip_profile=collections.defaultdict(int)
		
		#self.read_pos=[]
		#self.clip_count=[]
	
	def write(self,total_read):
		OUT=open(self._out_file1,'w')
		read_pos = []
		clip_count = []
		
		print >>OUT, "Position\tRead_Total\tRead_clipped"
		
		#max_pos = 0
		for i in self.soft_clip_profile:
			read_pos.append(str(i))
			clip_count.append(str(self.soft_clip_profile[i]))
			print >>OUT, str(i) + '\t' + str(total_read) + '\t' + str(self.soft_clip_profile[i])
		#	max_pos = i
		
		#max_pos += 1	
		if len(self.soft_clip_profile) > 0 :
			try :
				ROUT=open(self._out_file2,'w')
				#max_read_pos = max(read_pos)
				#sys.stderr.write(str(max_pos)+"\n")
				#sys.stderr.write(str(total_read)+"\n")
				print >>ROUT, "png(\"%s\",width=600,height=600,units='px')" % (self.plot_file)
				print >>ROUT, "read_pos=c(" + ','.join(read_pos) + ')'
				print >>ROUT, "count=c(" + ','.join(clip_count) + ')'
				#sys.stderr.write( 'plot(read_pos,1-(count/%d),xlim=c(0,%d),col="blue",main="Mappability Profile",xlab="Position of reads",ylab="Mappability",type="b")' % (total_read,max_read_pos))
				
				print >>ROUT, 'plot(read_pos,1-(count/'+str(total_read)+'),col="blue",xlab="Position of reads",ylab="Mappability",type="b")'
				print >>ROUT, "dev.state=dev.off()"
				ROUT.close()
			except :
				sys.stderr.write("Error in writing clipping plot script.\n")
				pass

		OUT.close()
	
	def set(self,cigar):
		cigar_str = bam_cigar.list2str(cigar)
		m = self.soft_p.findall(cigar_str)
		skip_part_of_read =1
		if len(m) > 0 :
			for j in m: 
				skip_part_of_read += sum([int(i) for i in self.read_part.findall(j[0])])
				for n in range(skip_part_of_read,(skip_part_of_read + int(j[1])+1)):
					self.soft_clip_profile[n]+=1
					#if n > 77 :
					#	sys.stderr.write(cigar_str+"\n")
					#	sys.exit(0)
				skip_part_of_read += int(j[1])
						
		
	
class ParseBAM:
	'''This class provides fuctions to parsing/processing/transforming SAM or BAM files. The input
	file could be either SAM or BAM format file'''
	
	multi_hit_tags=['H0','H1','H2','IH','NH']
	def __init__(self,inputFile):
		'''constructor. input could be bam or sam'''
		try:
			self.samfile = pysam.Samfile(inputFile,'rb')
			if len(self.samfile.header) ==0:
				print >>sys.stderr, "BAM/SAM file has no header section. Exit!"
				sys.exit(1)
			self.bam_format = True
		except:
			self.samfile = pysam.Samfile(inputFile,'r')
			if len(self.samfile.header) ==0:
				print >>sys.stderr, "BAM/SAM file has no header section. Exit!"
				sys.exit(1)
			self.bam_format = False

	def foundone(self,chrom,ranges, st, end):
		found = 0
		if chrom in ranges:
			found = len(ranges[chrom].find(st,end))
		return found

	def stat (self,exon_r,cds_exon,refbed,q_cut,outfile_fig,outfile,cds_exon_r, intron_r, utr_5_r, utr_3_r,intergenic_up_1kb_r,intergenic_down_1kb_r):
		
		res = Results.Results()
		'''Calculate mapping statistics'''
		intron_read = 0
		cds_exon_read = 0
		utr_5_read = 0
		utr_3_read = 0	
		intergenic_up1kb_read = 0
		intergenic_down1kb_read = 0
		unAssignFrags = 0
		intergenic_read = 0 
		
		prev_read_id = ""
		multi_read_flag = 0
		first_read_flag = 0
		second_read_flag = 0
		multi_reads = []
		multi_read1 = []
		multi_read2 = []
		
		clip_prof = Clipping_prof(outfile,outfile_fig)
		cov_prof = CoverageBody_prof(outfile,outfile_fig,exon_r)
		rDup_prof = ReadDup_prof(outfile,outfile_fig)
		inDist_prof = InnerDist_prof(outfile,outfile_fig)
		
		#bed_obj = BED.ParseBED(refbed)
 		ref_exons = []
 		
 		for exn in cds_exon:
 			ref_exons.append([exn[0].upper(), exn[1], exn[2]])		
		exon_bitsets = binned_bitsets_from_list(ref_exons)		
		
		try:
			while(1):
				
				aligned_read = self.samfile.next()
				#single end read
				cur_reads = []
				if not aligned_read.is_paired :					
					#count unmapped read
					if aligned_read.is_unmapped :
						res.unmapped_reads += 1
						res.total_reads += 1
						continue
									
					if aligned_read.qname == prev_read_id :
						multi_reads.append(aligned_read)
						continue
					else :				
						cur_read = None		
						if len(multi_reads) == 1 :
							cur_read = multi_reads[0]
							skip_read  = 0
							if cur_read.is_qcfail or cur_read.mapq < q_cut:			#skip QC fail read
								res.low_qual += 1
								skip_read = 1
							if cur_read.is_duplicate:		#skip duplicate read
								skip_read = 1	
													
							if skip_read == 1 :
								cur_read = None
							res.uniq_mapped_reads += 1
							res.total_reads += 1														
							
						if len(multi_reads) > 1:
							res.multi_mapped_reads += 1
							res.total_reads += 1
							cur_read = None
						
						if cur_read is not None :
							cur_reads.append(cur_read)
							if cur_read.is_reverse:
								res.reverse_read += 1
							else:
								res.forward_read += 1	
								
						if len(multi_reads) > 1:
							cur_reads.extend(multi_reads)
						multi_reads = []
						multi_reads.append(aligned_read)
						prev_read_id = aligned_read.qname
						
				#pair end read
				if aligned_read.is_paired :
					res.is_pairEnd = True					
					if aligned_read.qname == prev_read_id :
						if aligned_read.is_read1 :
							multi_read1.append(aligned_read)
						if aligned_read.is_read2 :
							multi_read2.append(aligned_read)						
					else :
						cur_read1 = None
						cur_read2 = None
						#uniq read
						if len(multi_read1) == 1 or len(multi_read2) == 1 :
							if len(multi_read1) == 1:
								cur_read1 = multi_read1[0]
							if len(multi_read2) == 1:
								cur_read2 = multi_read2[0]
							#unmapped read
							if (cur_read1 is None or cur_read1.is_unmapped) and (cur_read2 is None or cur_read2.is_unmapped) :
								res.unmapped_reads += 1
								
							if cur_read1 is None or cur_read1.is_unmapped :
								sys.stderr.write(cur_read1.qname+"\n")
								res.unmapped_read1 += 1
								cur_read1 = None
							else :
								res.mapped_read1 += 1
								if cur_read1.is_reverse :
									res.reverse_read += 1
								else :
									res.forward_read += 1
								if cur_read1.is_qcfail or cur_read1.mapq < q_cut:
									res.low_qual += 1
									cur_read1 = None
									
							if cur_read2 is None or cur_read2.is_unmapped :
								res.unmapped_read2 += 1
								sys.stderr.write(cur_read2.qname+"\n")
								cur_read2 = None
							else :
								res.mapped_read2 += 1
								if cur_read2.is_reverse :
									res.reverse_read += 1
								else :
									res.forward_read += 1
								if cur_read2.is_qcfail or cur_read2.mapq < q_cut:
									res.low_qual += 1
									cur_read2 = None
							
							#properly paired read
							if cur_read1 is not None and cur_read1.is_proper_pair :
								res.paired_reads += 1
								R_read1_ref = self.samfile.getrname(cur_read1.tid)
								R_read2_ref = self.samfile.getrname(cur_read1.rnext)
								if cur_read1.is_reverse and cur_read1.mate_is_reverse :
									res.mapped_minus_minus  += 1
								if cur_read1.is_reverse and not cur_read1.mate_is_reverse :
									res.mapped_minus_plus += 1
								if not cur_read1.is_reverse and cur_read1.mate_is_reverse :
									res.mapped_plus_minus += 1
								if not cur_read1.is_reverse and not cur_read1.mate_is_reverse :
									res.mapped_plus_plus += 1															
							
								if R_read1_ref != R_read2_ref:
									res.paired_diff_chrom += 1
									#cur_read1 = None
									#cur_read2 = None
								else :
									# computer inner dist
									inDist_prof.count(exon_bitsets, cur_read1, R_read1_ref.upper())								
								
							res.uniq_mapped_reads += 1
							res.total_reads += 1
							if cur_read1 is not None :
								cur_reads.append(cur_read1)
							if cur_read2 is not None :
								cur_reads.append(cur_read2)
								
						#multi-reads		
						if len(multi_read1) >1 or len(multi_read2) > 1 :
							
							res.multi_mapped_reads += 1
							res.total_reads += 1
							if len(multi_read1) > 1:
								res.mapped_read1 += 1
							if len(multi_read2) > 1 :
								res.mapped_read2 += 1						
							
						multi_read1 = []
						multi_read2 = []
						prev_read_id = aligned_read.qname
						if aligned_read.is_read1 :
							multi_read1.append(aligned_read)
						if aligned_read.is_read2 :
							multi_read2.append(aligned_read)
				
				if len(cur_reads) == 2 and res.is_pairEnd == True: #properly paired 
					cur_read = cur_reads[0]
					cur_read_mate = cur_reads[1]
					chrom1 = self.samfile.getrname(cur_read.tid)
					chrom1 = chrom1.upper()
					exons = bam_cigar.fetch_exon(chrom1, cur_read.pos, cur_read.cigar)
					
					chrom2 = self.samfile.getrname(cur_read_mate.tid)
					chrom2 = chrom2.upper()
					exons2 = bam_cigar.fetch_exon(chrom2, cur_read_mate.pos, cur_read_mate.cigar)				
					#calculate gene coverage				
					cov_prof.count2(exons,exons2,chrom1,chrom2)
					
				if len(cur_reads) == 1 : #either SE read or only one part is mapped of a PE
					cur_read = cur_reads[0]
					chrom = self.samfile.getrname(cur_read.tid)
					chrom = chrom.upper()
					exons = bam_cigar.fetch_exon(chrom, cur_read.pos, cur_read.cigar)
					#calculate gene coverage					
					cov_prof.count(exons,chrom)										
												
				for cur_read in cur_reads :
					clip_prof.set(cur_read.cigar)
					##########			
					chrom = self.samfile.getrname(cur_read.tid)
					chrom = chrom.upper()
					exons = bam_cigar.fetch_exon(chrom, cur_read.pos, cur_read.cigar)
					#total_Frags += len(exons)
				
					#calculate gene coverage
					#sys.stderr.write(cur_read.qname+"\n")					
					#cov_prof.count(exons,chrom,1.0/len(cur_reads))
					#--------------				
					#read duplicate rate
					RNA_read = cur_read.seq.upper()
					rDup_prof.count(RNA_read,chrom,cur_read.pos,exons)
					
					#---------------						
					for exn in exons:
						mid = int(exn[1]) + int((int(exn[2]) - int(exn[1])) / 2)
						if self.foundone(chrom, cds_exon_r, mid, mid) > 0:
								cds_exon_read += 1.0/len(cur_reads)
								continue
						elif self.foundone(chrom, utr_5_r, mid, mid) > 0 and self.foundone(chrom, utr_3_r, mid, mid) == 0:
								utr_5_read += 1.0/len(cur_reads)
								continue
						elif self.foundone(chrom, utr_3_r, mid, mid) > 0 and self.foundone(chrom, utr_5_r, mid, mid) == 0:
								utr_3_read += 1.0/len(cur_reads)
								continue
						elif self.foundone(chrom, utr_3_r, mid, mid) > 0 and self.foundone(chrom, utr_5_r, mid, mid) > 0:
								unAssignFrags += 1.0/len(cur_reads)
								continue
						elif self.foundone(chrom, intron_r, mid, mid) > 0:
								intron_read += 1.0/len(cur_reads)
								continue
						elif self.foundone(chrom, intergenic_up_1kb_r, mid, mid) > 0:
								intergenic_up1kb_read += 1.0/len(cur_reads)
								continue
						elif self.foundone(chrom, intergenic_down_1kb_r, mid, mid) > 0:
								intergenic_down1kb_read += 1.0/len(cur_reads)
								continue
						else:
							intergenic_read += 1.0/len(cur_reads)			
					
					introns = bam_cigar.fetch_intron('chr1', cur_read.pos, cur_read.cigar)
					if len(introns) == 0:
						res.noSplice += 1.0/len(cur_reads)
					if len(introns) >= 1:
						res.splice += 1.0/len(cur_reads)
								
		except StopIteration:
			pass
		
		res.splice = int(res.splice)
		res.noSplice = int(res.noSplice)
		
		prom_read = intergenic_up1kb_read  + intergenic_down1kb_read		  
		res.read_dist_plot_file = outfile_fig + ".read_distr.png"	
		  
		try :
			ROUT = open(outfile + '.read_distr.r', 'w')
			print >> ROUT, 'png("' + outfile_fig + '.read_distr.png",width=600,height=600,units="px")\n'
			print >> ROUT, "M=c(" + str(cds_exon_read) + "," + str(utr_5_read) + "," + str(utr_3_read) + "," + str(intron_read) + "," + str(prom_read) + "," + str(intergenic_read) + ")\n"
			
			print >> ROUT, "Mname=c('CDS','5UTR','3UTR','Intron','TSS_u/d1Kb','Others')\n"
			print >> ROUT, 'val = barplot(M,xlab="",space=1,ylab="Read Counts")\n'
			print >> ROUT, 'text(x=seq(val[1],val[6],by=2),y=rep(0,6),srt=60,adj=0,offset=2,pos=1,xpd=T,labels=Mname)\n'
			print >> ROUT, "dev.state = dev.off()"
			
			ROUT.close()
		except :
			sys.stderr.write("Error in writing plotting scripts.\n")
			pass
		
		clip_prof.write(res.total_reads)
		rDup_prof.write()
		inDist_prof.write()
		cov_prof.write(res.total_reads, refbed)
		
		res.insert_plot_file = inDist_prof.plot_file
		res.clipping_plot_file = clip_prof.plot_file
		res.read_dup_plot_file = rDup_prof.plot_file
		res.read_cov_plot_file = cov_prof.plot_file
		res.geneCount_file = cov_prof.outfile3
		
		return res

	
	def configure_experiment(self,refbed,sample_size = 200000, q_cut = 30):
		'''Given a BAM/SAM file, this function will try to guess the RNA-seq experiment:
			1) single-end or pair-end
			2) strand_specific or not
			3) if it is strand-specific, what's the strand_ness of the protocol
		'''
		
			#how many reads you want to sample
		count =0
		p_strandness=collections.defaultdict(int)
		s_strandness=collections.defaultdict(int)
		#load reference gene model
		gene_ranges={}
		print >>sys.stderr, "Reading reference gene model " + refbed + ' ...',
		for line in open(refbed,'r'):
			try:
				if line.startswith(('#','track','browser')):continue  
				# Parse fields from gene tabls
				fields = line.split()
				chrom	 = fields[0]
				tx_start  = int( fields[1] )
				tx_end	= int( fields[2] )
				geneName	  = fields[3]
				strand	= fields[5]
			except:
				print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line
				continue
			if chrom not in gene_ranges:
				gene_ranges[chrom]=Intersecter()
			gene_ranges[chrom].insert(tx_start,tx_end,strand)							
		print >>sys.stderr, "Done"		
		
		#read SAM/BAM file
		#current_pos = self.samfile.tell()
		print >>sys.stderr, "Loading SAM/BAM file ... ",
		try:
			while(1):
				if count >= sample_size:
					break
				aligned_read = self.samfile.next()
				if aligned_read.is_qcfail:			#skip low quanlity
					continue
				if aligned_read.is_duplicate:		#skip duplicate read
					continue
				if aligned_read.is_secondary:		#skip non primary hit
					continue
				if aligned_read.is_unmapped:		#skip unmap read
					continue		
				if aligned_read.mapq < q_cut:
					continue														
				
				chrom = self.samfile.getrname(aligned_read.tid)
				if aligned_read.is_paired:
					if aligned_read.is_read1:
						read_id = '1'
					if aligned_read.is_read2:
						read_id = '2'
					if aligned_read.is_reverse:
						map_strand = '-'
					else:
						map_strand = '+'
					readStart = aligned_read.pos
					readEnd = readStart + aligned_read.qlen
					if chrom in gene_ranges:
						if len(set(gene_ranges[chrom].find(readStart,readEnd)))==1:
							strand_from_gene = set(gene_ranges[chrom].find(readStart,readEnd)).pop()
							p_strandness[read_id + map_strand + strand_from_gene]+=1
							count += 1
				else:
					if aligned_read.is_reverse:
						map_strand = '-'
					else:
						map_strand = '+'					
					readStart = aligned_read.pos
					readEnd = readStart + aligned_read.qlen
					if chrom in gene_ranges:
						if len(set(gene_ranges[chrom].find(readStart,readEnd)))==1:
							strand_from_gene = set(gene_ranges[chrom].find(readStart,readEnd)).pop()
							s_strandness[map_strand + strand_from_gene]+=1
							count += 1

		except StopIteration:
			print >>sys.stderr, "Finished"		
		#self.samfile.seek(current_pos)
		
		print >>sys.stderr, "Total " + str(count) + " usable reads were sampled"
		protocol="unknown"
		strandness=None
		spec1=0.0
		spec2=0.0
		other=0.0
		if len(p_strandness) >0 and len(s_strandness) ==0 :
			protocol="PairEnd"
			#for k,v in p_strandness.items():
			#	print >>sys.stderr, k + '\t' + str(v)
			spec1= (p_strandness['1++'] + p_strandness['1--'] + p_strandness['2+-'] + p_strandness['2-+'])/float(sum(p_strandness.values()))
			spec2= (p_strandness['1+-'] + p_strandness['1-+'] + p_strandness['2++'] + p_strandness['2--'])/float(sum(p_strandness.values()))
			other = 1-spec1-spec2
			
		elif len(s_strandness) >0 and len(p_strandness) ==0 :
			protocol="SingleEnd"
			#for k,v in s_strandness.items():
			#	print  >>sys.stderr, k + '\t' + str(v)
			spec1 = (s_strandness['++'] + s_strandness['--'])/float(sum(s_strandness.values()))
			spec2 = (s_strandness['+-'] + s_strandness['-+'])/float(sum(s_strandness.values()))
			other = 1-spec1-spec2
		else:
			protocol="Mixture"
			spec1 = "NA"
			spec2 = "NA"
			other = "NA"
		return [protocol,spec1,spec2,other]

	def bamTowig(self,outfile,chrom_sizes, chrom_file,skip_multi=True,strand_rule=None,WigSumFactor=None,q_cut=30):
		"""Convert BAM/SAM file to wig file. chrom_size is dict with chrom as key and chrom_size as value
		strandRule should be determined from \"infer_experiment\". such as \"1++,1--,2+-,2-+\". When
		WigSumFactor is provided, output wig file will be normalized to this number """
		
		#strand_rule={'1+':'-','1-':'+','2+':'+','2-':'-'}
		strandRule={}
		if strand_rule is None:													# Not strand-specific
			pass																
		elif len(strand_rule.split(',')) ==4:									#PairEnd, strand-specific
			for i in strand_rule.split(','):strandRule[i[0]+i[1]]=i[2]
		elif len(strand_rule.split(',')) ==2:									#singeEnd, strand-specific
			for i in strand_rule.split(','):strandRule[i[0]]=i[1]
		else:
			print >>sys.stderr, "Unknown value of option :'strand_rule' " + strand_rule
			sys.exit(1)
		if len(strandRule) == 0:
			FWO = open(outfile + '.wig','w')
		else:
			FWO = open(outfile + '.Forward.wig','w')
			RVO = open(outfile + '.Reverse.wig','w')
		
		read_id=''
		
		for chr_name, chr_size in chrom_sizes.items():		#iterate each chrom
			try:
				self.samfile.fetch(chr_name,0,chr_size)
			except:
				print >>sys.stderr, "No alignments for " + chr_name + '. skipped'
				continue
			print >>sys.stderr, "Processing " + chr_name + " ..."
			if len(strandRule) == 0: FWO.write('variableStep chrom='+chr_name+'\n')
			else:
				FWO.write('variableStep chrom='+chr_name+'\n')
				RVO.write('variableStep chrom='+chr_name+'\n')
			Fwig = collections.defaultdict(int)
			Rwig = collections.defaultdict(int)
			alignedReads = self.samfile.fetch(chr_name,0,chr_size)		
			for aligned_read in alignedReads:
				if aligned_read.is_qcfail:continue			#skip low quanlity
				if aligned_read.is_duplicate:continue		#skip duplicate read
				if aligned_read.is_secondary:continue		#skip non primary hit
				if aligned_read.is_unmapped:continue		#skip unmap read
				
				if skip_multi:
					if aligned_read.mapq < q_cut:
						continue				
				if aligned_read.is_paired:
					if aligned_read.is_read1:read_id = '1'
					if aligned_read.is_read2:read_id = '2'
				
				if aligned_read.is_reverse:map_strand = '-'
				else:map_strand = '+'
				
				key = read_id + map_strand
				
				hit_st = aligned_read.pos
				for block in bam_cigar.fetch_exon(chr_name, hit_st, aligned_read.cigar): 
					for pos in range(block[1]+1,block[2]+1):	
						if len(strandRule) == 0: Fwig[pos] +=1.0	#this is NOT strand specific. everything into Fwig
						else:										#this is strand specific. separate Fwig and Rwig
							if strandRule[key] == '+':Fwig[pos] +=1.0
							if strandRule[key] == '-':Rwig[pos] -=1.0
			if WigSumFactor is None:	#not normalize	
				if len(strandRule) == 0:							#this is NOT strand specific.
					for pos in sorted (Fwig.keys()):
						print >>FWO, "%d\t%.2f" % (pos,Fwig[pos])
				else:
					for pos in sorted (Fwig.keys()):
						print >>FWO, "%d\t%.2f" % (pos,Fwig[pos])
					for pos in sorted (Rwig.keys()):
						print >>RVO, "%d\t%.2f" % (pos,Rwig[pos])
			else:						#normalize wig signal to WigSumFactor
				if len(strandRule) == 0:							#this is NOT strand specific.
					for pos in sorted (Fwig.keys()):
						print >>FWO, "%d\t%.2f" % (pos,Fwig[pos]*WigSumFactor)
				else:
					for pos in sorted (Fwig.keys()):
						print >>FWO, "%d\t%.2f" % (pos,Fwig[pos]*WigSumFactor)
					for pos in sorted (Rwig.keys()):
						print >>RVO, "%d\t%.2f" % (pos,Rwig[pos]*WigSumFactor)
		if len(strandRule) == 0:
			try:
				import subprocess
				print "Run " + "wigToBigWig " + outfile + '.wig ' + chrom_file + ' ' +  outfile + ".bw "
				subprocess.call("wigToBigWig " + outfile + '.wig ' + chrom_file + ' ' +  outfile + ".bw ",shell=True)
			except:
				print >>sys.stderr, "Cannot call \"wigToBigWig\". Conversion need to be done by yourself."
				pass
		else:
			try:
				import subprocess
				subprocess.call("wigToBigWig " + outfile + '.Forward.wig ' + chrom_file + ' ' +  outfile + ".Forward.bw ",shell=True)
				subprocess.call("wigToBigWig " + outfile + '.Reverse.wig ' + chrom_file + ' ' +  outfile + ".Reverse.bw ",shell=True)
			except:
				print >>sys.stderr, "Cannot call \"wigToBigWig\". Conversion need to be done by yourself."
				pass			

	def calWigSum(self,chrom_sizes, skip_multi=True):
		"""Calculate wigsum from BAM file"""
		
		print >>sys.stderr, "Calcualte wigsum ... "
		wigsum = 0.0
		read_id=''
		for chr_name, chr_size in chrom_sizes.items():		#iterate each chrom
			try:
				self.samfile.fetch(chr_name,0,chr_size)
			except:
				print >>sys.stderr, "No alignments for " + chr_name + '. skipped'
				continue
			print >>sys.stderr, "Processing " + chr_name + " ..."

			alignedReads = self.samfile.fetch(chr_name,0,chr_size)		
			for aligned_read in alignedReads:
				flag=0
				if aligned_read.is_qcfail:continue			#skip low quanlity
				if aligned_read.is_duplicate:continue		#skip duplicate read
				if aligned_read.is_secondary:continue		#skip non primary hit
				if aligned_read.is_unmapped:continue		#skip unmap read
				
				if skip_multi:
					if len(aligned_read.tags)>0:		#( ("NM", 1),("RG", "L1") )
						for i in aligned_read.tags:
							if i[0] in ParseBAM.multi_hit_tags and i[1] >1:
								flag=1						#multiple hit read
								break
					if flag==1:continue						#skip multiple map read		
				
				if aligned_read.is_paired:
					if aligned_read.is_read1:read_id = '1'
					if aligned_read.is_read2:read_id = '2'
				
				if aligned_read.is_reverse:map_strand = '-'
				else:map_strand = '+'
				
				key = read_id + map_strand
				
				hit_st = aligned_read.pos
				for block in bam_cigar.fetch_exon(chr_name, hit_st, aligned_read.cigar): 
					wigsum += (block[2] - block[1])
		return wigsum	

	def bam2fq(self,prefix, paired = True):
		"""Convert BAM into fastq files"""
		
		transtab = string.maketrans("ACGTNX","TGCANX")
		
		if paired:
			OUT1 = open(prefix + '.R1.fastq','w')
			OUT2 = open(prefix + '.R2.fastq','w')
			read1_count = 0
			read2_count = 0
		else:
			OUT = open(prefix + '.fastq','w')
			read_count = 0
		read_name = ''
		read_seq = ''
		read_qual = ''

		print >>sys.stderr, "Convert BAM file into fastq format ... ",
		try:
			while(1):
				aligned_read = self.samfile.next()
				read_name = aligned_read.qname
				read_seq = aligned_read.seq.upper()
				read_qual = aligned_read.qual
				if aligned_read.is_reverse:
					read_seq = read_seq.translate(transtab)[::-1]
					read_qual = read_qual[::-1]
				if paired:
					if aligned_read.is_read1:
						read1_count += 1
						if not read_name.endswith('/1'): 
							print >>OUT1, '@' + read_name + '/1'
						print >>OUT1, read_seq
						print >>OUT1, '+'
						print >>OUT1, read_qual
					if aligned_read.is_read2:
						read2_count += 1
						if not read_name.endswith('/2'): 
							print >>OUT2, '@' + read_name + '/2'
						print >>OUT2, read_seq
						print >>OUT2, '+'
						print >>OUT2, read_qual	
				else:		#single end
					read_count += 1
					print >>OUT, '@' + read_name
					print >>OUT, read_seq
					print >>OUT, '+'
					print >>OUT, read_qual	
					
		except StopIteration:
			print >>sys.stderr, "Done"
		if paired: 
			print >>sys.stderr, "read_1 count: %d" %  read1_count
			print >>sys.stderr, "read_2 count: %d" %  read2_count
		else:
			print >>sys.stderr, "read count: %d"  % read_count
					
	def calculate_rpkm(self,geneFile,outfile,strand_rule=None):
		'''calculate RPKM vaues. For single end RNA-seq, if it is strand specific, we assume that
		read plus mapped indicates a gene on plus strand.(similar to minus). 
		Advantages: works for both SAM and BAM
					works for both sorted and unsorted BAM/SAM file
					works for both index or unindexed BAM/SAM file
					much faster than indexing bam file
		Disadvantage: random access BAM file was disabled, thus large mount of RAM is required
		
		strand_rule: could be the following values:
			'1++,1--,2+-,2-+
			'1+-,1-+,2++,2--
			'++,--'
			'+-,-+'
			None
		'''
		
		strandRule={}
		if strand_rule is None:													# Not strand-specific
			pass																
		elif len(strand_rule.split(',')) ==4:									#PairEnd, strand-specific
			for i in strand_rule.split(','):strandRule[i[0]+i[1]]=i[2]
		elif len(strand_rule.split(',')) ==2:									#singeEnd, strand-specific
			for i in strand_rule.split(','):strandRule[i[0]]=i[1]
		else:
			print >>sys.stderr, "Unknown value of option :'strand_rule' " + strand_rule
			sys.exit(1)
		
		uniq_read=0
		total_tags=0
		plus_ranges={}
		minus_ranges={}
		unstrand_ranges={}
		
		rpkm_value={}
		
		RPKM_OUT = open(outfile,'w')
		if self.bam_format:print >>sys.stderr, "Load BAM file ... ",
		else:print >>sys.stderr, "Load SAM file ... ",
		
		#current_pos = self.samfile.tell()
		try:
			while(1):
				flag=0
				aligned_read = self.samfile.next()
				if aligned_read.is_qcfail:continue			#skip low quanlity					
				if aligned_read.is_duplicate:continue		#skip duplicate read
				if aligned_read.is_secondary:continue		#skip non primary hit
				if aligned_read.is_unmapped:continue		#skip unmap read

				if len(aligned_read.tags)>0:		#( ("NM", 1),("RG", "L1") )
					for i in aligned_read.tags:
						if i[0] in ParseBAM.multi_hit_tags and i[1] >1:
							flag=1						#multiple hit read
							break
				if flag==1:continue						#skip multiple map read		
				
				uniq_read +=1
				
				if aligned_read.is_paired:
					if aligned_read.is_read1:read_id = '1'
					if aligned_read.is_read2:read_id = '2'
				else:
					read_id = ''
				if aligned_read.is_reverse:map_strand = '-'
				else:map_strand = '+'
				
				strand_key = read_id + map_strand
				
				chrom = self.samfile.getrname(aligned_read.tid).upper()
				hit_st = aligned_read.pos
				exon_blocks = bam_cigar.fetch_exon(chrom, hit_st, aligned_read.cigar)
				total_tags += len(exon_blocks)
						
				#construct bitset
				if strand_rule is not None:	
					if strandRule[strand_key] == '+':
						for block in exon_blocks:
							mid = block[1] + int((block[2] - block[1])/2)
							if chrom not in plus_ranges:plus_ranges[chrom] = Intersecter()
							plus_ranges[chrom].add_interval( Interval( mid,mid+1 ) )
					elif strandRule[strand_key] == '-':
						for block in exon_blocks:
							mid = block[1] + int((block[2] - block[1])/2)
							if chrom not in minus_ranges:minus_ranges[chrom] = Intersecter()	
							minus_ranges[chrom].add_interval( Interval( mid,mid+1 ) )
				elif strand_rule is None:	
					for block in exon_blocks:
						mid = block[1] + int((block[2] - block[1])/2)
						if chrom not in unstrand_ranges:unstrand_ranges[chrom] = Intersecter()
						unstrand_ranges[chrom].add_interval( Interval( mid,mid+1 ) )
					
		except StopIteration:
			print >>sys.stderr, "Done"
		#self.samfile.seek(current_pos)
		print >>RPKM_OUT, "#Total uniquely mapped reads = " + str(uniq_read)
		print >>RPKM_OUT, "#Total fragments = " + str(total_tags)
		print >>sys.stderr, "Assign reads to "+ geneFile + '...',
		for line in open(geneFile,'r'):
			try:
				if line.startswith('#'):continue
				if line.startswith('track'):continue
				if line.startswith('browser'):continue   
				# Parse fields from gene tabls
				fields = line.split()
				chrom	 = fields[0].upper()
				tx_start  = int( fields[1] )
				tx_end	= int( fields[2] )
				geneName	  = fields[3]
				strand	= fields[5].replace(" ","_")
				
				exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
				exon_starts = map((lambda x: x + tx_start ), exon_starts)
				exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
				exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends)
				exon_sizes = map(int,fields[10].rstrip(',\n').split(','))
				intron_starts = exon_ends[:-1]
				intron_ends=exon_starts[1:]
				key='\t'.join((chrom.lower(),str(tx_start),str(tx_end),geneName,'0',strand))
			except:
				print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line,
				continue

					
			mRNA_count=0
			mRNA_len=sum(exon_sizes)
				
			if (strand_rule is not None) and (strand == '-'):
				intronNum=len(intron_starts)
				exonNum=len(exon_starts)
				
				# assign reads to intron	
				for st,end in zip(intron_starts,intron_ends):
					if chrom in minus_ranges:
						hits= len(minus_ranges[chrom].find(st,end))
						RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(total_tags))) +'\n')
						intronNum -= 1
				# assign reads to exon				
				for st,end in zip(exon_starts,exon_ends):
					if chrom in minus_ranges:
						hits= len(minus_ranges[chrom].find(st,end))					
						RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(total_tags))) +'\n')
						exonNum -= 1
						mRNA_count += hits
				try:
					RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(mRNA_count) + "\t" + strand + '\t' +  str(mRNA_count*1000000000.0/(mRNA_len*total_tags)) +'\n')
				except:
					RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(0) + "\t" + strand + '\t' +  str(0) +'\n')
			elif (strand_rule is not None) and (strand == '+'):
				intronNum=1
				exonNum=1
				for st,end in zip(intron_starts,intron_ends):
					if chrom in plus_ranges:
						hits= len(plus_ranges[chrom].find(st,end))
						RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(total_tags))) +'\n')
						intronNum += 1	
				for st,end in zip(exon_starts,exon_ends):
					if chrom in plus_ranges:
						hits= len(plus_ranges[chrom].find(st,end))
						RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(total_tags))) +'\n')
						exonNum += 1		
						mRNA_count += hits
				try:
					RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(mRNA_count) + "\t" + strand + '\t' +  str(mRNA_count*1000000000.0/(mRNA_len*total_tags)) +'\n')
				except:
					RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(0) + "\t" + strand + '\t' +  str(0) +'\n')
			elif strand_rule is None:
				intronNum=1
				exonNum=1
				for st,end in zip(intron_starts,intron_ends):
					if chrom in unstrand_ranges:
						hits= len(unstrand_ranges[chrom].find(st,end))
						RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(total_tags))) +'\n')
						intronNum += 1	
				for st,end in zip(exon_starts,exon_ends):
					if chrom in unstrand_ranges:
						hits= len(unstrand_ranges[chrom].find(st,end))
						RPKM_OUT.write(chrom.lower() + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_exon_" + str(exonNum) + "\t" + str(hits) + "\t" + strand + '\t' +  str(hits*1000000000.0/((end-st)*(total_tags))) +'\n')
						exonNum += 1		
						mRNA_count += hits
				try:
					RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(mRNA_count) + "\t" + strand + '\t' +  str(mRNA_count*1000000000.0/(mRNA_len*total_tags)) +'\n')
				except:
					RPKM_OUT.write(chrom.lower() + "\t" + str(tx_start) + "\t" + str(tx_end) + "\t" + geneName + "_mRNA" + "\t" + str(0) + "\t" + strand + '\t' +  str(0) +'\n')
		print >>sys.stderr, "Done"

	def readsNVC(self,outfile=None,nx=True, q_cut = 30):
		'''for each read, calculate nucleotide frequency vs position'''
		if outfile is None:
			outfile1 = self.fileName + ".NVC.xls"
			outfile2 = self.fileName +".NVC_plot.r"
		else:
			outfile1 = outfile + ".NVC.xls"
			outfile2 = outfile +".NVC_plot.r"
		FO=open(outfile1,'w')
		RS=open(outfile2,'w')
		PPcount=0
		
		transtab = string.maketrans("ACGTNX","TGCANX")
		base_freq=collections.defaultdict(int)
		a_count=[]
		c_count=[]
		g_count=[]
		t_count=[]
		n_count=[]
		x_count=[]
		if self.bam_format:print >>sys.stderr, "Read BAM file ... ",
		else:print >>sys.stderr, "Read SAM file ... ",

		try:
			while(1):
				aligned_read = self.samfile.next()
				if aligned_read.mapq < q_cut: continue
				#if aligned_read.is_unmapped:continue	#skip unmapped read
				#if aligned_read.is_qcfail:continue	#skip low quality

				RNA_read = aligned_read.seq.upper()		
				if aligned_read.is_reverse:
					RNA_read = RNA_read.translate(transtab)[::-1]
				for i,j in enumerate(RNA_read):
					key = str(i) + j
					base_freq[key] += 1
		except StopIteration:
			print >>sys.stderr, "Done"
		
		print >>sys.stderr, "generating data matrix ..."
		print >>FO, "Position\tA\tC\tG\tT\tN\tX"
		for i in xrange(len(RNA_read)):
			print  >>FO, str(i) + '\t',
			print  >>FO, str(base_freq[str(i) + "A"]) + '\t',
			a_count.append(str(base_freq[str(i) + "A"]))
			print  >>FO, str(base_freq[str(i) + "C"]) + '\t',
			c_count.append(str(base_freq[str(i) + "C"]))
			print  >>FO, str(base_freq[str(i) + "G"]) + '\t',
			g_count.append(str(base_freq[str(i) + "G"]))
			print  >>FO, str(base_freq[str(i) + "T"]) + '\t',
			t_count.append(str(base_freq[str(i) + "T"]))
			print  >>FO, str(base_freq[str(i) + "N"]) + '\t',
			n_count.append(str(base_freq[str(i) + "N"]))
			print  >>FO, str(base_freq[str(i) + "X"]) + '\t'
			x_count.append(str(base_freq[str(i) + "X"]))
		FO.close()
		
		#generating R scripts
		print >>sys.stderr, "generating R script  ..."
		print >>RS, "position=c(" + ','.join([str(i) for i in xrange(len(RNA_read))]) + ')'
		print >>RS, "A_count=c(" + ','.join(a_count) + ')'
		print >>RS, "C_count=c(" + ','.join(c_count) + ')'
		print >>RS, "G_count=c(" + ','.join(g_count) + ')'
		print >>RS, "T_count=c(" + ','.join(t_count) + ')'
		print >>RS, "N_count=c(" + ','.join(n_count) + ')'
		print >>RS, "X_count=c(" + ','.join(x_count) + ')'
		
		if nx:
			print >>RS, "total= A_count + C_count + G_count + T_count + N_count + X_count"
			print >>RS, "ym=max(A_count/total,C_count/total,G_count/total,T_count/total,N_count/total,X_count/total) + 0.05"
			print >>RS, "yn=min(A_count/total,C_count/total,G_count/total,T_count/total,N_count/total,X_count/total)"
			
			print >>RS, 'png(\"%s\")' % (outfile +".NVC_plot.png")
			print >>RS, 'plot(position,A_count/total,type="o",pch=20,ylim=c(yn,ym),col="dark green",xlab="Position of Read",ylab="Nucleotide Frequency")'
			print >>RS, 'lines(position,T_count/total,type="o",pch=20,col="red")'
			print >>RS, 'lines(position,G_count/total,type="o",pch=20,col="blue")'
			print >>RS, 'lines(position,C_count/total,type="o",pch=20,col="cyan")'
			print >>RS, 'lines(position,N_count/total,type="o",pch=20,col="black")'		
			print >>RS, 'lines(position,X_count/total,type="o",pch=20,col="grey")'	
			print >>RS, 'legend('+ str(len(RNA_read)-10) + ',ym,legend=c("A","T","G","C","N","X"),col=c("dark green","red","blue","cyan","black","grey"),lwd=2,pch=20,text.col=c("dark green","red","blue","cyan","black","grey"))'
			print >>RS, "dev.state = dev.off()"
		else:
			print >>RS, "total= A_count + C_count + G_count + T_count"
			print >>RS, "ym=max(A_count/total,C_count/total,G_count/total,T_count/total) + 0.05"
			print >>RS, "yn=min(A_count/total,C_count/total,G_count/total,T_count/total)"
		
			print >>RS, 'png(\"%s\")' % (outfile +".NVC_plot.png")
			print >>RS, 'plot(position,A_count/total,type="o",pch=20,ylim=c(yn,ym),col="dark green",xlab="Position of Read",ylab="Nucleotide Frequency")'
			print >>RS, 'lines(position,T_count/total,type="o",pch=20,col="red")'
			print >>RS, 'lines(position,G_count/total,type="o",pch=20,col="blue")'
			print >>RS, 'lines(position,C_count/total,type="o",pch=20,col="cyan")'
			print >>RS, 'legend('+ str(len(RNA_read)-10) + ',ym,legend=c("A","T","G","C"),col=c("dark green","red","blue","cyan"),lwd=2,pch=20,text.col=c("dark green","red","blue","cyan"))'
			print >>RS, "dev.state = dev.off()"
		
		RS.close()
		#self.f.seek(0)

	def readsQual_boxplot(self,outfile,shrink=1000, q_cut=30):
		'''calculate phred quality score for each base in read (5->3)'''

		output = outfile + ".qual.r"
		FO=open(output,'w')

		if self.bam_format:print >>sys.stderr, "Read BAM file ... ",
		else:print >>sys.stderr, "Read SAM file ... ",

		quality = collections.defaultdict(dict)	#read_pos=>quality score=>count
		q_max = -1
		q_min = 10000
		q_list=[]
		i_box={}	#key is read postion,value is 
		try:
			while(1):
				aligned_read = self.samfile.next()
				if aligned_read.mapq < q_cut: continue
				#if aligned_read.is_unmapped:continue	#skip unmapped read
				#if aligned_read.is_qcfail:continue		#skip low quality
				
				qual_str = aligned_read.qqual
				read_len = aligned_read.rlen
				if aligned_read.is_reverse:
					qual_str = qual_str[::-1]

				for i,j in enumerate(qual_str):
					q=ord(j)-33
					if q > q_max: q_max = q
					if q < q_min: q_min = q
					try:
						quality[i][q] += 1
					except:
						quality[i][q] = 1
		except StopIteration:
			print >>sys.stderr, "Done"
		
		for p in range(0,read_len):
			#print str(p) + ':',
			val=[]
			occurrence=[]
			for q in range(q_min,q_max+1):
				if quality.has_key(p) and quality[p].has_key(q):
					val.append(str(q))				
					occurrence.append(str(quality[p][q]))	
					q_list.append(str(quality[p][q]))
				else:
					q_list.append(str(0))
			i_box[p] = 'rep(c(' + ','.join(val) + '),times=c(' + ','.join(occurrence) + ')/' + str(shrink)+ ')'
		
		
		#generate R script for boxplot
		print >>FO, "png(\'%s\')" % (outfile + ".qual.boxplot.png")
		for i in sorted(i_box):
			print >>FO,'p'+str(i) + '<-' + i_box[i]
		print >>FO, 'boxplot(' + ','.join(['p'+str(i) for i in i_box]) + ',xlab=\"Position of Read(5\'->3\')\",ylab=\"Phred Quality Score\",outline=F' + ')'
		print >>FO,"dev.state = dev.off()"
		
		
		#generate R script for heatmap
		print >>FO, '\n'
		print >>FO, "png(\'%s\')" % (outfile + ".qual.heatmap.png")
		print >>FO, "qual=c(" + ','.join(q_list)  + ')'
		print >>FO, "mat=matrix(qual,ncol=%s,byrow=F)" % (read_len)
		print >>FO, 'Lab.palette <- colorRampPalette(c("blue", "orange", "red3","red2","red1","red"), space = "rgb",interpolate=c(\'spline\'))'
		print >>FO, "heatmap(mat,Rowv=NA,Colv=NA,xlab=\"Position of Read\",ylab=\"Phred Quality Score\",labRow=seq(from=%s,to=%s),col = Lab.palette(256),scale=\"none\" )" % (q_min,q_max)
		print >>FO, 'dev.state = dev.off()'
		
	def readGC(self,outfile=None, q_cut=30):
		'''GC content distribution of reads'''
		if outfile is None:
			outfile1 = self.fileName + ".GC.xls"
			outfile2 = self.fileName +".GC_plot.r"
		else:
			outfile1 = outfile + ".GC.xls"
			outfile2 = outfile + ".GC_plot.r"
		FO=open(outfile1,'w')
		RS=open(outfile2,'w')
		
		gc_hist=collections.defaultdict(int)	#key is GC percent, value is count of reads

		if self.bam_format:print >>sys.stderr, "Read BAM file ... ",
		else:print >>sys.stderr, "Read SAM file ... ",

		try:
			while(1):
				aligned_read = self.samfile.next()
				if aligned_read.is_unmapped:continue	#skip unmapped read
				if aligned_read.is_qcfail:continue		#skip low quality
				if aligned_read.mapq < q_cut: continue
				RNA_read = aligned_read.seq.upper()		
				gc_percent = "%4.2f" % ((RNA_read.count('C') + RNA_read.count('G'))/(len(RNA_read)+0.0)*100)
				#print gc_percent
				gc_hist[gc_percent] += 1
		except StopIteration:
			print >>sys.stderr, "Done"
		
		print >>sys.stderr, "writing GC content ..."	
		print >>FO, "GC%\tread_count"
		for i in gc_hist.keys():
			print >>FO, i + '\t' + str(gc_hist[i])
			
		print >>sys.stderr, "writing R script ..."
		print >>RS, "png(\"%s\")" % (outfile +  ".GC_plot.png")
		print >>RS, 'gc=rep(c(' + ','.join([i for i in gc_hist.keys()]) + '),' + 'times=c(' + ','.join([str(i) for i in gc_hist.values()]) + '))'
		print >>RS, 'hist(gc,probability=T,breaks=%d,xlab="GC content (%%)",ylab="Density of Reads",border="blue",main="")' % 100
		#print >>RS, "lines(density(gc),col='red')"
		print >>RS ,"dev.state = dev.off()"		
		#self.f.seek(0)
		

	def annotate_junction(self,refgene,outfile,min_intron=50, q_cut=30):
		'''Annotate splicing junctions in BAM or SAM file. Note that a (long) read might have multiple splicing
		events  (splice multiple times), and the same splicing events can be consolidated into a single
		junction'''
		
		out_file = outfile + ".junction.xls"
		out_file2 = outfile + ".junction_plot.r"
		if refgene is None:
			print >>sys.stderr, "You must provide reference gene model in bed format."
			sys.exit(1)
		OUT = open(out_file,'w')
		ROUT = open(out_file2,'w')
		
		#reading reference gene model
		refIntronStarts=collections.defaultdict(dict)
		refIntronEnds=collections.defaultdict(dict)	
		total_junc =0
		novel35_junc =0
		novel3or5_junc =0
		known_junc =0
		splicing_events=collections.defaultdict(int)	
		
		print >>sys.stderr, "Reading reference bed file: ",refgene, " ... ",
		for line in open(refgene,'r'):
			if line.startswith(('#','track','browser')):continue  
		   	# Parse fields from gene tabls
			fields = line.split()
			if(len(fields)<12):
				print >>sys.stderr, "Invalid bed line (skipped):",line,
				continue
			chrom	 = fields[0].upper()
			tx_start = int( fields[1] )
			tx_end   = int( fields[2] )
			if int(fields[9] ==1):
				continue		
			
			exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
			exon_starts = map((lambda x: x + tx_start ), exon_starts)
			exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
			exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends);   
			intron_start = exon_ends[:-1]
			intron_end=exon_starts[1:]
			for i_st,i_end in zip (intron_start, intron_end):
				refIntronStarts[chrom][i_st] =i_st
				refIntronEnds[chrom][i_end] =i_end			
		print >>sys.stderr,"Done"
		
		#reading input SAM file
		if self.bam_format:print >>sys.stderr, "Load BAM file ... ",
		else:print >>sys.stderr, "Load SAM file ... ",

		try:
			while(1):
				aligned_read = self.samfile.next()
				if aligned_read.is_qcfail:continue			#skip low quanlity					
				if aligned_read.is_duplicate:continue		#skip duplicate read
				if aligned_read.is_secondary:continue		#skip non primary hit
				if aligned_read.is_unmapped:continue		#skip unmap read
				if aligned_read.mapq < q_cut:continue

				chrom = self.samfile.getrname(aligned_read.tid).upper()
				hit_st = aligned_read.pos
				intron_blocks = bam_cigar.fetch_intron(chrom, hit_st, aligned_read.cigar)			
				if len(intron_blocks)==0:
					continue
				for intrn in intron_blocks:
					total_junc +=1
					if intrn[2] - intrn[1] < min_intron:continue
					splicing_events[intrn[0] + ":" + str(intrn[1]) + ":" + str(intrn[2])] += 1
					if (refIntronStarts[chrom].has_key(intrn[1]) and refIntronEnds[chrom].has_key(intrn[2])):
						known_junc +=1																		#known both
					elif (not refIntronStarts[chrom].has_key(intrn[1]) and not refIntronEnds[chrom].has_key(intrn[2])):
						novel35_junc +=1																
					else:
						novel3or5_junc +=1
		except StopIteration:
			print >>sys.stderr, "Done"
		
		print "total = " + str(total_junc)
		if total_junc == 0:
			print >>sys.stderr, "No splice junction found."
			sys.exit()
		#self.f.seek(0)
		
		print >>ROUT, 'png(\"%s\")' % (outfile + ".splice_events.png")
		print >>ROUT, "events=c(" + ','.join([str(i*100.0/total_junc) for i in (novel3or5_junc,novel35_junc,known_junc)])+ ')'
		print >>ROUT, 'pie(events,col=c(2,3,4),init.angle=30,angle=c(60,120,150),density=c(70,70,70),main="splicing events",labels=c("partial_novel %d%%","complete_novel %d%%","known %d%%"))' % (round(novel3or5_junc*100.0/total_junc),round(novel35_junc*100.0/total_junc),round(known_junc*100.0/total_junc))
		print >>ROUT, "dev.state = dev.off()"
		
		print >>sys.stderr, "\n==================================================================="
		print >>sys.stderr, "Total splicing  Events:\t" + str(total_junc)
		print >>sys.stderr, "Known Splicing Events:\t" + str(known_junc)
		print >>sys.stderr, "Partial Novel Splicing Events:\t" + str(novel3or5_junc)
		print >>sys.stderr, "Novel Splicing Events:\t" + str(novel35_junc)
		
		#reset variables
		total_junc =0
		novel35_junc =0
		novel3or5_junc =0
		known_junc =0
		
		print >>OUT, "chrom\tintron_st(0-based)\tintron_end(1-based)\tread_count\tannotation"
		for i in splicing_events:
			total_junc += 1
			(chrom, i_st, i_end) = i.split(":")
			print >>OUT, '\t'.join([chrom.replace("CHR","chr"),i_st,i_end]) + '\t' + str(splicing_events[i]) + '\t',
			i_st = int(i_st)
			i_end = int(i_end)
			if (refIntronStarts[chrom].has_key(i_st) and refIntronEnds[chrom].has_key(i_end)):
				print >>OUT, "annotated"
				known_junc +=1
			elif (not refIntronStarts[chrom].has_key(i_st) and not refIntronEnds[chrom].has_key(i_end)):
				print >>OUT, 'complete_novel'
				novel35_junc +=1
			else:
				print >>OUT, 'partial_novel'
				novel3or5_junc +=1
		
		if total_junc ==0:
			print >>sys.stderr, "No splice read found"
			sys.exit(1)
		print >>sys.stderr, "\nTotal splicing  Junctions:\t" + str(total_junc)
		print >>sys.stderr, "Known Splicing Junctions:\t" + str(known_junc)
		print >>sys.stderr, "Partial Novel Splicing Junctions:\t" + str(novel3or5_junc)
		print >>sys.stderr, "Novel Splicing Junctions:\t" + str(novel35_junc)
		print >>sys.stderr, "\n==================================================================="
		
		print >>ROUT, 'png(\"%s\")' % (outfile + ".splice_junction.png")
		print >>ROUT, "junction=c(" + ','.join([str(i*100.0/total_junc) for i in (novel3or5_junc,novel35_junc,known_junc,)])+ ')'
		print >>ROUT, 'pie(junction,col=c(2,3,4),init.angle=30,angle=c(60,120,150),density=c(70,70,70),main="splicing junctions",labels=c("partial_novel %d%%","complete_novel %d%%","known %d%%"))' % (round(novel3or5_junc*100.0/total_junc),round(novel35_junc*100.0/total_junc),round(known_junc*100.0/total_junc))
		print >>ROUT, "dev.state = dev.off()"
		#print >>ROUT, "mat=matrix(c(events,junction),byrow=T,ncol=3)"
		#print >>ROUT, 'barplot(mat,beside=T,ylim=c(0,100),names=c("known","partial\nnovel","complete\nnovel"),legend.text=c("splicing events","splicing junction"),ylab="Percent")'

	def saturation_junction(self,refgene,outfile=None,sample_start=5,sample_step=5,sample_end=100,min_intron=50,recur=1, q_cut=30):
		'''check if an RNA-seq experiment is saturated in terms of detecting known splicing junction'''
		
		out_file = outfile + ".junctionSaturation_plot.r"
		if refgene is None:
			print >>sys.stderr, "You must provide reference gene model in bed format."
			sys.exit(1)
		
		OUT = open(out_file,'w')


		#reading reference gene 
		knownSpliceSites= set()
		chrom_list=set()
		print >>sys.stderr, "reading reference bed file: ",refgene, " ... ",
		for line in open(refgene,'r'):
			if line.startswith(('#','track','browser')):continue  
			fields = line.split()
			if(len(fields)<12):
				print >>sys.stderr, "Invalid bed line (skipped):",line,
				continue
			chrom	 = fields[0].upper()
			chrom_list.add(chrom)
			tx_start = int( fields[1] )
			tx_end   = int( fields[2] )
			if int(fields[9] ==1):
				continue		
			
			exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
			exon_starts = map((lambda x: x + tx_start ), exon_starts)
			exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
			exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends);   
			intron_start = exon_ends[:-1]
			intron_end=exon_starts[1:]
			for st,end in zip (intron_start, intron_end):
				knownSpliceSites.add(chrom + ":" + str(st) + "-" + str(end))
		print >>sys.stderr,"Done! Total "+str(len(knownSpliceSites)) + " known splicing junctions."


		#read SAM file
		samSpliceSites=[]
		intron_start=[]
		intron_end=[]
		uniqSpliceSites=collections.defaultdict(int)

		if self.bam_format:print >>sys.stderr, "Load BAM file ... ",
		else:print >>sys.stderr, "Load SAM file ... ",
		try:
			while(1):
				aligned_read = self.samfile.next()
				chrom = self.samfile.getrname(aligned_read.tid).upper()
				if chrom not in chrom_list:
					continue				
				if aligned_read.is_qcfail:continue			#skip low quanlity					
				if aligned_read.is_duplicate:continue		#skip duplicate read
				if aligned_read.is_secondary:continue		#skip non primary hit
				if aligned_read.is_unmapped:continue		#skip unmap read
				if aligned_read.mapq < q_cut: continue
				
				hit_st = aligned_read.pos
				intron_blocks = bam_cigar.fetch_intron(chrom, hit_st, aligned_read.cigar)			
				if len(intron_blocks)==0:
					continue
				for intrn in intron_blocks:
					if intrn[2] - intrn[1] < min_intron:continue
					samSpliceSites.append(intrn[0] + ":" + str(intrn[1]) + "-" + str(intrn[2]))
		except StopIteration:
			print >>sys.stderr, "Done"
		
		print >>sys.stderr, "shuffling alignments ...",
		random.shuffle(samSpliceSites)
		print >>sys.stderr, "Done"
				
		#resampling
		SR_num = len(samSpliceSites)
		sample_size=0
		all_junctionNum = 0	
		known_junc=[]
		all_junc=[]
		unknown_junc=[]
		#=========================sampling uniquely mapped reads from population
		tmp=range(sample_start,sample_end,sample_step)
		tmp.append(100)
		for pertl in tmp:	#[5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95,100]
			knownSpliceSites_num = 0
			index_st = int(SR_num * ((pertl - sample_step)/100.0))
			index_end = int(SR_num * (pertl/100.0))
			if index_st < 0: index_st = 0
			sample_size += index_end -index_st
			
			print >>sys.stderr, "sampling " + str(pertl) +"% (" + str(sample_size) + ") splicing reads.",
			
			#all splice juntion
			for i in range(index_st, index_end):
				uniqSpliceSites[samSpliceSites[i]] +=1	
			all_junctionNum = len(uniqSpliceSites.keys())
			all_junc.append(str(all_junctionNum))
			print >>sys.stderr, str(all_junctionNum) + " splicing junctions.",
			
			#known splice junction
			known_junctionNum = 0
			for sj in uniqSpliceSites:
				if sj in knownSpliceSites and uniqSpliceSites[sj] >= recur:
					known_junctionNum +=1
			print >>sys.stderr, str(known_junctionNum) + " known splicing junctions.",
			known_junc.append(str(known_junctionNum))
			
			#unknown splice junction
			unknown_junctionNum = 0
			for sj in uniqSpliceSites:
				if sj not in knownSpliceSites:
					unknown_junctionNum +=1
			unknown_junc.append(str(unknown_junctionNum))
			print >>sys.stderr, str(unknown_junctionNum) + " novel splicing junctions."
			
		#for j in uniq_SJ:
			#print >>OUT, j + "\t" + str(uniq_SJ[j])
		print >>OUT, "png(\'%s\')" % (outfile + '.junctionSaturation_plot.png')
		print >>OUT, "x=c(" + ','.join([str(i) for i in tmp]) + ')'
		print >>OUT, "y=c(" + ','.join(known_junc) + ')'
		print >>OUT, "z=c(" + ','.join(all_junc) + ')'
		print >>OUT, "w=c(" + ','.join(unknown_junc) + ')'
		print >>OUT, "m=max(%d,%d,%d)" % (int(int(known_junc[-1])/1000), int(int(all_junc[-1])/1000),int(int(unknown_junc[-1])/1000))
		print >>OUT, "n=min(%d,%d,%d)" % (int(int(known_junc[0])/1000), int(int(all_junc[0])/1000),int(int(unknown_junc[0])/1000))
		print >>OUT, "plot(x,z/1000,xlab='percent of total reads',ylab='Number of splicing junctions (x1000)',type='o',col='blue',ylim=c(n,m))"
		print >>OUT, "points(x,y/1000,type='o',col='red')"
		print >>OUT, "points(x,w/1000,type='o',col='green')"
		print >>OUT, 'legend(5,%d, legend=c("All junctions","known junctions", "novel junctions"),col=c("blue","red","green"),lwd=1,pch=1)' % int(int(all_junc[-1])/1000)
		print >>OUT, "dev.state = dev.off()"

	def saturation_RPKM(self,refbed,outfile,sample_start=5,sample_step=5,sample_end=100,skip_multi=True, strand_rule=None, q_cut=30):
		'''for each gene, check if its RPKM (epxresion level) has already been saturated or not'''
		
		if refbed is None:
			print >>sys.stderr,"You must specify a bed file representing gene model\n"
			exit(0)
		rpkm_file = outfile + ".eRPKM.xls"
		raw_file = outfile + ".rawCount.xls"
		
		RPKM_OUT = open(rpkm_file,'w')
		RAW_OUT = open(raw_file ,'w')
		
		ranges={}
		totalReads=0
		cUR_num = 0	#number of fragements
		cUR_plus = 0
		cUR_minus = 0
		block_list_plus = []	#non-spliced read AS IS, splicing reads were counted multiple times
		block_list_minus = []
		block_list = []
		strandRule = {}
				
		if strand_rule is None:													# Not strand-specific
			pass																
		elif len(strand_rule.split(',')) ==4:									#PairEnd, strand-specific
			for i in strand_rule.split(','):strandRule[i[0]+i[1]]=i[2]
		elif len(strand_rule.split(',')) ==2:									#singeEnd, strand-specific
			for i in strand_rule.split(','):strandRule[i[0]]=i[1]
		else:
			print >>sys.stderr, "Unknown value of: 'strand_rule' " +  strand_rule
			sys.exit(1)	


		#read SAM or BAM
		if self.bam_format:print >>sys.stderr, "Load BAM file ... ",
		else:print >>sys.stderr, "Load SAM file ... ",
		try:
			while(1):
				aligned_read = self.samfile.next()
				if aligned_read.is_qcfail:continue			#skip low quanlity					
				if aligned_read.is_duplicate:continue		#skip duplicate read
				if aligned_read.is_secondary:continue		#skip non primary hit
				if aligned_read.is_unmapped:continue		#skip unmap read
				
				if skip_multi:
					if aligned_read.mapq < q_cut:
						continue			
				chrom = self.samfile.getrname(aligned_read.tid).upper()
				
				#determine read_id and read_strand
				if aligned_read.is_paired:						#pair end
					if aligned_read.is_read1:read_id = '1'
					if aligned_read.is_read2:read_id = '2'
				else:read_id = ''								#single end
			
				if aligned_read.is_reverse:map_strand = '-'
				else:map_strand = '+'				
				strand_key = read_id + map_strand				#used to determine if a read should assign to gene(+) or gene(-)

				hit_st = aligned_read.pos
				exon_blocks = bam_cigar.fetch_exon(chrom, hit_st, aligned_read.cigar)	
				cUR_num += len(exon_blocks)
				
				#strand specific
				if strand_rule is not None:
					if strandRule[strand_key] == '+': cUR_plus += len(exon_blocks)
					if strandRule[strand_key] == '-': cUR_minus += len(exon_blocks)
					for exn in exon_blocks:
						if strandRule[strand_key] == '+': block_list_plus.append(exn[0] + ":" + str(exn[1] + (exn[2]-exn[1])/2 ))
						if strandRule[strand_key] == '-': block_list_minus.append(exn[0] + ":" + str(exn[1] + (exn[2]-exn[1])/2 ))
				#Not strand specific
				else:			
					for exn in exon_blocks:
						block_list.append(exn[0] + ":" + str(exn[1] + (exn[2]-exn[1])/2 ))
		except StopIteration:
			print >>sys.stderr, "Done"
		
		
		print >>sys.stderr, "shuffling alignments ...",
		random.shuffle(block_list_plus)
		random.shuffle(block_list_minus)
		random.shuffle(block_list)
		print >>sys.stderr, "Done"
		
		
		ranges_plus={}
		ranges_minus={}
		ranges={}
		sample_size=0
		RPKM_table=collections.defaultdict(list)
		rawCount_table=collections.defaultdict(list)
		RPKM_head=['#chr','start','end','name','score','strand']

		tmp=range(sample_start,sample_end,sample_step)
		tmp.append(100)
		#=========================sampling uniquely mapped reads from population
		for pertl in tmp:	#[5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95,100]
			percent_st = (pertl-sample_step)/100.0
			percent_end = pertl/100.0
			if percent_st < 0: percent_st = 0
			sample_size = cUR_num * percent_end
			RPKM_head.append(str(pertl) + '%')
			
			if strand_rule is not None:
				print >>sys.stderr, "sampling " + str(pertl) +"% (" + str(int(cUR_plus * percent_end)) + ") forward strand fragments ..."
				for i in block_list_plus[int(cUR_plus*percent_st):int(cUR_plus*percent_end)]:
					(chr,coord) = i.split(':')
					if chr not in ranges_plus:ranges_plus[chr] = Intersecter()
					ranges_plus[chr].add_interval( Interval( int(coord), int(coord)+1 ) )								
				
				print >>sys.stderr, "sampling " + str(pertl) +"% (" + str(int(cUR_minus * percent_end)) + ") reverse strand fragments ..."			
				for i in block_list_minus[int(cUR_minus*percent_st):int(cUR_minus*percent_end)]:
					(chr,coord) = i.split(':')
					if chr not in ranges_minus:ranges_minus[chr] = Intersecter()				
					ranges_minus[chr].add_interval( Interval( int(coord), int(coord)+1 ) )						
			
			else:
				print >>sys.stderr, "sampling " + str(pertl) +"% (" + str(int(sample_size)) + ") fragments ..."
				for i in block_list[int(cUR_num*percent_st):int(cUR_num*percent_end)]:
					(chr,coord) = i.split(':')
					if chr not in ranges:ranges[chr] = Intersecter()						
					ranges[chr].add_interval( Interval( int(coord), int(coord)+1 ) )														

			#========================= calculating RPKM based on sub-population
			print >>sys.stderr, "assign reads to transcripts in " + refbed + ' ...'
			for line in open(refbed,'r'):
				try:
					if line.startswith(('#','track','browser')):continue  
					# Parse fields from gene tabls
					fields = line.split()
					chrom	 = fields[0].upper()
					tx_start  = int( fields[1] )
					tx_end	= int( fields[2] )
					geneName	  = fields[3]
					strand	= fields[5]
					exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
					exon_starts = map((lambda x: x + tx_start ), exon_starts)
					exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
					exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends)
					exon_sizes = map(int,fields[10].rstrip(',\n').split(','))
					key='\t'.join((chrom.lower(),str(tx_start),str(tx_end),geneName,'0',strand))
				except:
					print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line
					continue
				mRNA_count=0	#we need to initializ it to 0 for each gene
				mRNA_len=sum(exon_sizes)
				for st,end in zip(exon_starts,exon_ends):
					#if chrom in ranges:
					if strand_rule is not None:
						if (strand == '+') and (chrom in ranges_plus): mRNA_count += len(ranges_plus[chrom].find(st,end))	
						if (strand == '-') and (chrom in ranges_minus): mRNA_count += len(ranges_minus[chrom].find(st,end))
					else:
						if chrom in ranges:
							mRNA_count += len(ranges[chrom].find(st,end))
				if mRNA_len ==0:
					print >>sys.stderr, geneName + " has 0 nucleotides. Exit!"
					sys.exit(1)
				if sample_size == 0:
					print >>sys.stderr, "Too few reads to sample. Exit!"
					sys.exit(1)
				mRNA_RPKM = (mRNA_count * 1000000000.0)/(mRNA_len * sample_size)
				RPKM_table[key].append(str(mRNA_RPKM))
				rawCount_table[key].append(str(mRNA_count))
			print >>sys.stderr, ""

		#self.f.seek(0)
		print >>RPKM_OUT, '\t'.join(RPKM_head)
		print >>RAW_OUT, '\t'.join(RPKM_head)
		for key in RPKM_table:
			print >>RPKM_OUT, key + '\t',
			print >>RPKM_OUT, '\t'.join(RPKM_table[key])
			print >>RAW_OUT, key + '\t',
			print >>RAW_OUT, '\t'.join(rawCount_table[key])		


	def shuffle_RPKM(self,refbed,outfile,sample_percentage=0.5,shuffle_times=50,skip_multi=True, strand_rule=None):
		'''for each gene, check if its RPKM (epxresion level) has already been saturated or not'''
		
		if refbed is None:
			print >>sys.stderr,"You must specify a bed file representing gene model\n"
			exit(0)
		rpkm_file = outfile + ".eRPKM.xls"
		raw_file = outfile + ".rawCount.xls"
		
		RPKM_OUT = open(rpkm_file,'w')
		RAW_OUT = open(raw_file ,'w')
		
		ranges={}
		totalReads=0
		cUR_num = 0	#number of fragements
		cUR_plus = 0
		cUR_minus = 0
		block_list_plus = []	#non-spliced read AS IS, splicing reads were counted multiple times
		block_list_minus = []
		block_list = []
		strandRule = {}
				
		if strand_rule is None:													# Not strand-specific
			pass																
		elif len(strand_rule.split(',')) ==4:									#PairEnd, strand-specific
			for i in strand_rule.split(','):strandRule[i[0]+i[1]]=i[2]
		elif len(strand_rule.split(',')) ==2:									#singeEnd, strand-specific
			for i in strand_rule.split(','):strandRule[i[0]]=i[1]
		else:
			print >>sys.stderr, "Unknown value of: 'strand_rule' " +  strand_rule
			sys.exit(1)	


		#read SAM or BAM
		if self.bam_format:print >>sys.stderr, "Load BAM file ... ",
		else:print >>sys.stderr, "Load SAM file ... ",
		try:
			while(1):
				flag=0
				aligned_read = self.samfile.next()
				if aligned_read.is_qcfail:continue			#skip low quanlity					
				if aligned_read.is_duplicate:continue		#skip duplicate read
				if aligned_read.is_secondary:continue		#skip non primary hit
				if aligned_read.is_unmapped:continue		#skip unmap read
				
				if skip_multi:
					if len(aligned_read.tags)>0:		#( ("NM", 1),("RG", "L1") )
						for i in aligned_read.tags:
							if i[0] in ParseBAM.multi_hit_tags and i[1] >1:
								flag=1						#multiple hit read
								break
				if flag==1:continue						#skip multiple map read		
				
				chrom = self.samfile.getrname(aligned_read.tid).upper()
				
				#determine read_id and read_strand
				if aligned_read.is_paired:						#pair end
					if aligned_read.is_read1:read_id = '1'
					if aligned_read.is_read2:read_id = '2'
				else:read_id = ''								#single end
			
				if aligned_read.is_reverse:map_strand = '-'
				else:map_strand = '+'				
				strand_key = read_id + map_strand				#used to determine if a read should assign to gene(+) or gene(-)

				hit_st = aligned_read.pos
				exon_blocks = bam_cigar.fetch_exon(chrom, hit_st, aligned_read.cigar)	
				cUR_num += len(exon_blocks)
				
				#strand specific
				if strand_rule is not None:
					if strandRule[strand_key] == '+': cUR_plus += len(exon_blocks)
					if strandRule[strand_key] == '-': cUR_minus += len(exon_blocks)
					for exn in exon_blocks:
						if strandRule[strand_key] == '+': block_list_plus.append(exn[0] + ":" + str(exn[1] + (exn[2]-exn[1])/2 ))
						if strandRule[strand_key] == '-': block_list_minus.append(exn[0] + ":" + str(exn[1] + (exn[2]-exn[1])/2 ))
				#Not strand specific
				else:			
					for exn in exon_blocks:
						block_list.append(exn[0] + ":" + str(exn[1] + (exn[2]-exn[1])/2 ))
		except StopIteration:
			print >>sys.stderr, "Done"
				
		RPKM_table=collections.defaultdict(list)
		rawCount_table=collections.defaultdict(list)
		RPKM_head=['#chr','start','end','name','score','strand']
		
		iter_times=0
		#=========================sampling uniquely mapped reads from population
		for x in range(0,shuffle_times+1):
			print >>sys.stderr, "Shuffle " + str(iter_times) + " times"
			iter_times += 1
			if iter_times == shuffle_times:
				sample_percent = 1
			else:
				sample_percent = sample_percentage
			ranges_plus={}
			ranges_minus={}
			ranges={}
			if strand_rule is not None:
				for i in random.sample(block_list_plus, int(cUR_plus * sample_percent)):
					(chr,coord) = i.split(':')
					if chr not in ranges_plus:ranges_plus[chr] = Intersecter()
					ranges_plus[chr].add_interval( Interval( int(coord), int(coord)+1 ) )								
				
				for i in random.sample(block_list_minus, int(cUR_minus * sample_percent)):
					(chr,coord) = i.split(':')
					if chr not in ranges_minus:ranges_minus[chr] = Intersecter()				
					ranges_minus[chr].add_interval( Interval( int(coord), int(coord)+1 ) )						
			
			else:
				for i in random.sample(block_list,int(cUR_num * sample_percent)):
					(chr,coord) = i.split(':')
					if chr not in ranges:ranges[chr] = Intersecter()						
					ranges[chr].add_interval( Interval( int(coord), int(coord)+1 ) )														

			#========================= calculating RPKM based on sub-population
			print >>sys.stderr, "assign reads to transcripts in " + refbed + ' ...'
			for line in open(refbed,'r'):
				try:
					if line.startswith(('#','track','browser')):continue  
					# Parse fields from gene tabls
					fields = line.split()
					chrom	 = fields[0].upper()
					tx_start  = int( fields[1] )
					tx_end	= int( fields[2] )
					geneName	  = fields[3]
					strand	= fields[5]
					exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
					exon_starts = map((lambda x: x + tx_start ), exon_starts)
					exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
					exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends)
					exon_sizes = map(int,fields[10].rstrip(',\n').split(','))
					key='\t'.join((chrom.lower(),str(tx_start),str(tx_end),geneName,'0',strand))
				except:
					print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line
					continue
				mRNA_count=0	#we need to initializ it to 0 for each gene
				mRNA_len=sum(exon_sizes)
				for st,end in zip(exon_starts,exon_ends):
					#if chrom in ranges:
					if strand_rule is not None:
						if (strand == '+') and (chrom in ranges_plus): mRNA_count += len(ranges_plus[chrom].find(st,end))	
						if (strand == '-') and (chrom in ranges_minus): mRNA_count += len(ranges_minus[chrom].find(st,end))
					else:
						if chrom in ranges:
							mRNA_count += len(ranges[chrom].find(st,end))
				if mRNA_len ==0:
					print >>sys.stderr, geneName + " has 0 nucleotides. Exit!"
					sys.exit(1)
				if cUR_num * sample_percentage == 0:
					print >>sys.stderr, "Too few reads to sample. Exit!"
					sys.exit(1)
				mRNA_RPKM = (mRNA_count * 1000000000.0)/(mRNA_len * (cUR_num * sample_percentage))
				RPKM_table[key].append(str(mRNA_RPKM))
				rawCount_table[key].append(str(mRNA_count))
			print >>sys.stderr, ""

		#self.f.seek(0)
		print >>RPKM_OUT, '\t'.join(RPKM_head)
		print >>RAW_OUT, '\t'.join(RPKM_head)
		for key in RPKM_table:
			print >>RPKM_OUT, key + '\t',
			print >>RPKM_OUT, '\t'.join(RPKM_table[key])
			print >>RAW_OUT, key + '\t',
			print >>RAW_OUT, '\t'.join(rawCount_table[key])		
		
	def fetchAlignments(self,chr,st,end):
		'''fetch alignment from sorted BAM file based on chr, st, end
		Note: BAM file must be indexed'''
		try:
			a=self.samfile.fetch(chr,st,end)
			return a
		except:
			return None
		

def print_bits_as_bed( bits ):
	end = 0
	while 1:
		start = bits.next_set( end )
		if start == bits.size: break
		end = bits.next_clear( start )
		print "%d\t%d" % ( start, end )
