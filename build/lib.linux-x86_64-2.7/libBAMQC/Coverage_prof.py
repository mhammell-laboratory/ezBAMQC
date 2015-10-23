#!/usr/bin/env python
'''manipulate BAM/SAM file.'''
'''Created by Ying Jin Aug. 18, 2015 '''

#import built-in modules
import os,sys, re, string
import warnings

import collections
import math
import sets
import random


def percentile_list(N):
        """
            Find the percentile of a list of values.
            @parameter N - is a list of values. Note N MUST BE already sorted.
            @return - the list of percentile of the values
            """
        if not N: return None
        per_list=[]
        for i in range(0,101):
            k = (len(N)-1) * i/100.0
            f = math.floor(k)
            c = math.ceil(k)
            if f == c:
                    per_list.append( int(N[int(k)])  )
            else:
                    d0 = N[int(f)] * (c-k)
                    d1 = N[int(c)] * (k-f)
                    per_list.append(int(round(d0+d1)))
        return per_list



class CoverageBody_prof:
    def __init__(self,outfile_data,outfile_fig,geneIdx):
        self.outfile1 = outfile_data + ".geneBodyCoverage_plot.r"
        self.outfile2 = outfile_data + ".geneBodyCoverage.txt"
        self.plot_file = outfile_fig + ".geneBodyCoverage.png"
        self.plot_file2 = outfile_fig + ".TransCoverage.png"
        self.outfile4 = outfile_data + ".TransCoverage.r"
        self.outfile3 = outfile_data + ".geneAbundance.txt"
        
        self.ranges = {}
        self.frag_num = 0
        self.geneCounts_list = [0] * len(geneIdx.getFeatures())
        self.gene_percentile_base = {}
        genes = geneIdx.getFeatures()

        self.geneToidx = dict(zip(genes,range(0,len(genes))))
        
        #for g in self.geneToidx :
        #    sys.stderr.write(g + "\t" + str(self.geneToidx[g]) + "\n")
                
        for g in genes :
            gene_all_base = []
            for st_end_str in geneIdx.get_exon_pair(g) :
                #sys.stderr.write(st_end_str + "\n")
                (st,end,strand) = st_end_str.split(':')
                gene_all_base.extend(range(int(st),int(end)))
                if strand == "-" :
                    gene_all_base.sort(reverse = True)
                
            mRNA_len = len(gene_all_base)
            if mRNA_len <100:
                    continue
            perc_list = percentile_list(gene_all_base)
            self.gene_percentile_base[g] = (perc_list,[0]*len(perc_list))


    
    def write(self,totalReads):
        coverage=collections.defaultdict(int)
        total_exons = 0
        zero_exons = 0

        for g in self.gene_percentile_base:
            (idx,counts) = self.gene_percentile_base[g]
            
            for i in range(len(counts)) :
                total_exons += 1
                coverage[i] += counts[i]
                if counts[i] == 0 :
                    zero_exons += 1
                    
        x_coord=[]
        y_coord=[]
        try :
            OUT2 = open(self.outfile2,'w')
            OUT1 = open(self.outfile1,'w')
            OUT3 = open(self.outfile3,'w')
            OUT4 = open(self.outfile4,'w')
                
            print >>OUT4, "png(\'%s\',width=500,height=500,units='px')" % (self.plot_file2)
            print >>OUT4, "a=c("+ str(self.geneCounts_list).strip('[]') +")"
            print >>OUT4, "Fn = ecdf(a)"
            print >>OUT4, "max_x = round(log(max(knots(Fn)),2),0)"
            print >>OUT4, "xx = c(0,2^seq(0,max_x,by=2))"
            print >>OUT4, "y=Fn(xx)"
            print >>OUT4, 'plot(x=xx,y=y,type="b",col="blue",pch=20,xlab="Number of Reads",ylab="Cumulative proportion of Genes")\n'
            print >>OUT4, "dev.state = dev.off()"
            OUT4.close()
            
            print >>OUT2, "Total reads: " + str(totalReads)
            print >>OUT2, "Fragment number: " + str(self.frag_num)
            print >>OUT2, "percentile\tcount"
            
            for i in range(len(coverage)):
                x_coord.append(str(i))
                y_coord.append(str(coverage[i]))
                print >>OUT2, str(i) + '\t' + str(coverage[i])
            
            print >>OUT1, "png(\'%s\',width=500,height=500,units='px')" % (self.plot_file)
            print >>OUT1, "x=c("+','.join(x_coord)+')'
            print >>OUT1, "y=c(" + ','.join(y_coord) + ')'
            print >>OUT1, "smoothsp = smooth.spline(x,y,spar=0.35)"
            print >>OUT1, 'plot(smoothsp,type="l",col="blue",xlab="Percentile of Gene Body (5\'->3\')",ylab="Number of read",xlim=c(0,100))\n'
            print >>OUT1, "dev.state = dev.off()"
            
            OUT1.close()
            OUT2.close()
            
            if len(self.geneCounts_list) >0 :
                #sys.stderr.write(str(len(self.geneCounts_list))+"\n")
                for key in sorted(self.geneToidx.keys()) :
                    #sys.stderr.write(key+"\n")
                    print >>OUT3, key +"\t"+str(self.geneCounts_list[self.geneToidx[key]])
                
            OUT3.close()
            
        except :
            sys.stderr.write("Error in writing gene body coverage data.\n")
            pass
        return (total_exons,zero_exons)

#gene abundance
    def count(self,gene,exon_blocks,exon_blocks2):
        if gene != "" :
            self.frag_num += len(exon_blocks)+len(exon_blocks2)
            self.geneCounts_list[self.geneToidx[gene]] += 1
            self.__per_base_count(gene,exon_blocks,exon_blocks2)


#update gene percentile counts
    def __per_base_count(self,gene,exon_blocks1,exon_blocks2) :

        ovp_base = []
        
        exons1_sorted = sorted(exon_blocks1,key=lambda tup:tup[0])
        for (st,end) in exons1_sorted :
            ovp_base.extend(range(st,end))
        
        exons2_sorted = sorted(exon_blocks2,key=lambda tup:tup[0])
        for (st,end) in exons2_sorted :
            ovp_base.extend(range(st,end))
        
        if gene in self.gene_percentile_base :
            (gene_perc_base,gene_perc_base_cnt) = self.gene_percentile_base[gene]
            for i in range(len(gene_perc_base)) :
                if gene_perc_base[i] in ovp_base :
                    gene_perc_base_cnt[i] += 1





