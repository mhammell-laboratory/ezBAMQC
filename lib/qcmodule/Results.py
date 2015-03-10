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
#changes to the paths

class Results :
    
    def __init__(self):
        self.filename = ""
        self.is_pairEnd = False
        self.clipping_plot_file = ""
        self.read_cov_plot_file = ""
        self.insert_plot_file = ""
        self.read_dist_plot_file = ""
        self.read_dup_plot_file = ""
        self.geneCount_file = ""
        
        self.no_clipping = False
        self.total_reads = 0
        self.uniq_mapped_reads = 0
        self.multi_mapped_reads = 0
        self.unmapped_reads = 0
        self.low_qual = 0
        self.pcr_dup = 0
        
        self.unmapped_read1 = 0
        self.unmapped_read2 = 0
        self.mapped_read1 = 0
        self.mapped_read2 = 0 
        self.forward_read = 0
        self.reverse_read = 0
        self.paired_reads = 0
        
        self.mapped_plus_minus = 0
        self.mapped_plus_plus = 0
        self.mapped_minus_plus = 0
        self.mapped_minus_minus = 0
        
        self.ins_read = 0
        self.del_read = 0
        
        self.noSplice = 0
        self.splice = 0
        self.paired_diff_chrom = 0
       
        
'''    def stat2 (self):
        
        results = Results()
        #Calculate mapping statistics        
        total_read=0
        paired_reads = 0
        
        pcr_duplicate =0
        low_qual =0
        secondary_hit =0
        
        unmapped_read1=0
        mapped_read1=0
        reverse_read1=0
        forward_read1=0
        
        unmapped_read2=0
        mapped_read2=0
        reverse_read2=0
        forward_read2=0
        
        _numSplitHit =0
        _numMonoHit =0
        _numInsertion =0
        _numDeletion =0

        minus_minus=0
        minus_plus =0
        plus_minus=0
        plus_plus=0
        paired=True
        
        unmap_SE=0
        map_SE=0
        reverse_SE=0
        forward_SE=0
        
        multi_reads = 0
        
        prev_read_id = ""
        prev_read1_id = ""
        prev_read2_id = ""
        multi_read_flag = 0
        
        for line in self.f:
            line=line.rstrip()
            if line.startswith('@'):continue                #skip head lines    
            if ParseSAM._reExpr2.match(line):continue        #skip blank lines
            
        #    total_read +=1
            field = line.split()
            flagCode=string.atoi(field[1])
            
            read_id = field[0]
            
            # paired end reads  
            if (flagCode & 0x0001 !=0):                                #This is paired end sequencing                
                if (flagCode & 0x0002 != 0 ) :#both pair properly mapped
                    # read in the next line (part of current pair)
                    line2 = self.f.readline()
                    line2 = line2.rstrip()
                    field2 = line2.split()
                    
                    flagCode2 = string.atoi(field2[1])
                    if read_id == prev_read1_id :
                        if multi_read_flag == 0 :
                            multi_read_flag = 1
                            multi_read += 1
                            
                    if read_id != prev_read1_id :
                        total_read += 1
                        paired_reads += 1
                        if  (flagCode & 0x0400 !=0):                            #PCR or optical duplicate
                            pcr_duplicate +=1
                            continue
                        if  (flagCode & 0x0200 !=0):                            #Low quality
                            low_qual +=1
                            continue 
                        mapped_read1 += 1
                        mapped_read2 += 1
                        multi_read_flag = 0
                        prev_read1_id = read_id
                        prev_read2_id = read_id
                        if (flagCode & 0x0010 != 0):reverse_read1 +=1
                        if (flagCode & 0x0010 == 0):forward_read1 +=1
                        if (flagCode2 & 0x0010 != 0):reverse_read2 +=1
                        if (flagCode2 & 0x0010 == 0):forward_read2 +=1                          
                
                else : # only one is mapped or both are not mapped    
                       if  (flagCode & 0x0400 !=0):                            #PCR or optical duplicate
                              total_read += 1
                              pcr_duplicate +=1
                              continue
                       if  (flagCode & 0x0200 !=0):                            #Low quality
                              low_qual +=1
                              total_read += 1
                              continue
                       if (flagCode & 0x0040 != 0):                        #1st read
                           if read_id != prev_read1_id :
                               if read_id != prev_read2_id :
                                   total_read += 1
                                   multi_read_flag = 0
                            if (flagCode & 0x0004 != 0):unmapped_read1 +=1
                            if (flagCode & 0x0004 == 0):
                                  mapped_read1 +=1
                                  prev_read1_id = read_id
                                  if (flagCode & 0x0010 != 0):reverse_read1 +=1
                                  if (flagCode & 0x0010 == 0):forward_read1 +=1
                           else :
                                  if multi_read_flag ==0 :
                                      multi_read_flag = 1
                                      multi_reads += 1
                           if (flagCode & 0x0080 != 0):                        #2nd read
                            if read_id != prev_read2_id :
                                if read_id != prev_read1_id :
                                    total_read += 1
                                    multi_read_flag = 0
                                if (flagCode & 0x0004 != 0):unmapped_read2 +=1
                                if (flagCode & 0x0004 == 0):
                                    mapped_read2 +=1
                                    prev_read2_id = read_id
                                
                                   if (flagCode & 0x0010 != 0):reverse_read2 +=1
                                   if (flagCode & 0x0010 == 0):forward_read2 +=1
                            else :
                                   if multi_read_flag ==0 :
                                       multi_read_flag = 1
                                       multi_reads += 1                                
#            if     (flagCode & 0x0010 != 0 and flagCode & 0x0020 != 0):
#                    minus_minus +=1
#                if     (flagCode & 0x0010 != 0 and flagCode & 0x0020 == 0):
#                    minus_plus +=1
#                if     (flagCode & 0x0010 == 0 and flagCode & 0x0020 != 0):
#                    plus_minus +=1
#                if     (flagCode & 0x0010 == 0 and flagCode & 0x0020 == 0):
#                    plus_plus +=1        
#                    if  (flagCode & 0x0100 !=0):                            #Not primary alignment
#                        secondary_hit +=1
#                        continue
#                    if (len(ParseSAM._splicedHit_pat.findall(field[5]))>1):_numSplitHit +=1            #Splicing mapped reads                            
#                    if (len(ParseSAM._splicedHit_pat.findall(field[5]))==1):_numMonoHit +=1            #mono mapped reads                    
#                    if (ParseSAM._insertionHit_pat.search(field[5])):_numInsertion +=1                #insertion in reads
#                    if (ParseSAM._deletionHit_pat.search(field[5])):_numDeletion +=1                #deletion in reads
#*****************This is single end sequencing***********************
           if (flagCode & 0x0001 ==0):                                
                        paired=False
                        if  (flagCode & 0x0400 !=0):                            #PCR or optical duplicate
                                    total_read += 1
                                    pcr_duplicate +=1
                                    continue
                        if  (flagCode & 0x0200 !=0):                            #Low quality
                                low_qual +=1
                                total_read += 1
                                continue
                        
                        if (flagCode & 0x0004 != 0):
                            unmap_SE +=1
                            
                        if (flagCode & 0x0004 == 0):
                            if read_id != prev_read_id :
                                map_SE +=1
                                total_read += 1

                            else :
                                if multi_read_flag ==0 :
                                    multi_read_flag = 1
                                    multi_read += 1'''
                                    
                        #if (flagCode & 0x0010 != 0):
                        #    reverse_SE +=1
                        #if (flagCode & 0x0010 == 0):
                        #    forward_SE +=1
                
'''        if paired:        
                results.total_reads = total_read 
                results.mapped_read1 = mapped_read1
                results.mapped_read2 = mapped_read2
                results.uniq_mapped_reads = total_read - low_qual - pcr_dup - multi_reads 
                results.multi_mapped_reads = multi_reads 
                results.low_qual = low_qual
                results.pcr_dup = pcr_duplicate
                
                results.unmapped_read1 = unmapped_read1
                results.unmapped_read2 = unmapped_read2

                results.forward_read1 = forward_read1
                results.reverse_read1 = reverse_read1
                results.forward_read2 = forward_read2
                results.reverse_read2 = reverse_read2
                results.paired_reads = paired_reads
                
        else:
                results.total_reads = total_read 
                results.uniq_mapped_reads = map_SE - multi_reads
                results.multi_mapped_reads = multi_reads 
                results.low_qual = low_qual
                results.pcr_dup = pcr_duplicate'''
        
            
#        return results
 
