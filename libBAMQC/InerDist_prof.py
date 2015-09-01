#!/usr/bin/env python
'''manipulate BAM/SAM file.'''

#import built-in modules
import os,sys
import re
import string
import collections
import math


class InnerDist_prof:
    def __init__(self,outfile_data,outfile_fig,sample_size=100000,low_bound=0,up_bound=1000,step=10):
        
        #self.out_file1 = outfile_data + ".inner_distance.txt"
        self.out_file2 = outfile_data + ".inner_distance_freq.txt"
        self.out_file3 = outfile_data + ".inner_distance_plot.r"
        self.plot_file = outfile_fig+".inner_distance_plot.png"
        
        self.pair_num = 0
        self.sample_size = sample_size
        self.step = step
        self.window_left_bound = range(low_bound,up_bound,step)
        self.counts = [0] *(up_bound-low_bound+1)
        
    def write(self):
        '''estimate the inner distance of mRNA pair end fragment. fragment size = insert_size + 2 x read_length'''
        if self.pair_num == 0 :
            return
            
        try :
            #FO=open(self.out_file1,'w')
            FQ=open(self.out_file2,'w')
            RS=open(self.out_file3,'w')
        except:
            sys.stderr.write("Error in creating file inner distance frequency files. \n")
            pass
        
        sizes = []
        counts = []
        
        for st in self.window_left_bound:
            sizes.append(str(st + self.step/2))
            count = sum(self.counts[st:(st+self.step)])
            counts.append(str(count))
            print >>FQ, str(st) + '\t' + str(st+self.step-1) +'\t' + str(count)
        
        #plot_file = outfile_fig+".inner_distance_plot.png"
        print >>RS, "png(\"%s\",width=500,height=500,units='px')" % (self.plot_file)

        print >>RS, 'fragsize=rep(c(' + ','.join(sizes) + '),' + 'times=c(' + ','.join(counts) + '))'
        print >>RS, 'frag_sd = round(sd(fragsize),digits=0)'
        print >>RS, 'frag_mean = round(mean(fragsize),digits=0)'
        print >>RS, 'frag_median = round(median(fragsize),digits=0)'
        print >>RS, 'write(c("Mean insert size",frag_mean), stdout())'
        print >>RS, 'write(c("Median insert size",frag_median), stdout())'
        print >>RS, 'write(c("Standard deviation",frag_sd), stdout())'
        print >>RS, 'hist(fragsize,probability=T,breaks=%d,xlab="mRNA insert size (bp)",main=paste(c("Mean=",frag_mean,";","SD=",frag_sd),collapse=""),border="blue")' % len(self.window_left_bound)
        print >>RS, "lines(density(fragsize,bw=%d),col='red')" % (2*self.step)
        print >>RS ,"dev.state = dev.off()"
        #FO.close()
        FQ.close()
        RS.close()
        
        
        #    def count(self,exon_bitsets,aligned_read,chrom): cds_exon
    def count(self,geneIdx,gene,aligned_read,chrom,intron_blocks,strand):

        if self.pair_num > self.sample_size:
            return
        self.pair_num +=1

        splice_intron_size=0            

        read1_len = aligned_read.infer_query_length()
            
        read1_start = aligned_read.pos
        read2_start = aligned_read.mpos


        for intron in intron_blocks:
            splice_intron_size += intron[2] - intron[1]
                
        read1_end = read1_start + read1_len + splice_intron_size        

        inner_distance = read2_start - read1_end +1
        
        if inner_distance > 0 and inner_distance < len(self.counts):
            exons = []
            if gene is not None :
                exons = geneIdx.find_exons(chrom,read1_end,read2_start,strand,gene)
            
            if len(exons) > 0 :
                exon_starts =  map(lambda tup: tup[0],exons)
                exon_ends =  map(lambda tup: tup[1],exons)
                size = sum(map(lambda x, y : x - y + 1, exon_ends,exon_starts))
                
                if size <= inner_distance and size > 1:
                        self.counts[size] += 1
                else:
                    self.counts[inner_distance] += 1
                        
            else:
                    self.counts[inner_distance] += 1
            

            
