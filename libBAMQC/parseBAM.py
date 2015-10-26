#!/usr/bin/env python
'''manipulate BAM/SAM file.'''

#import built-in modules
import os,sys, time
import re
import string
import warnings
import collections
import math
import sets
import random
import pysam

from libBAMQC.Results import *
from libBAMQC.Mappability import *
from libBAMQC.Coverage_prof import *
from libBAMQC.InerDist_prof import *
from libBAMQC.ReadDup_prof import *


def fetch_intron( st, cigar,format):
    ''' fetch intron regions defined by cigar. st must be zero based
        return list of tuple of (st, end)
        '''
    #match = re.compile(r'(\d+)(\D)')
    chrom_st = st
    if format == "BAM" :
        chrom_st += 1
    
    intron_bound =[]

    for (c,s) in cigar:    #code and size
        
        if c==0:        #match
            chrom_st += s
        elif c==1:        #insertion to ref
            continue
        elif c==2:        #deletion to ref
            chrom_st += s
        elif c==3:        #gap or intron
            intron_bound.append([chrom_st,chrom_st+s])
            chrom_st += s
        elif c==4:        #soft clipping. We do NOT include soft clip as part of exon
            chrom_st += s
        else:
            continue
        
    return intron_bound

def fetch_exon(chrom, st, cigar,format):
    ''' fetch exon regions defined by cigar. st must be zero based
        return list of tuple of (st, end)
        '''
    #match = re.compile(r'(\d+)(\D)')
    chrom_st = st
    if format == "BAM" :
        chrom_st += 1
    exon_bound =[]

    for c,s in cigar:    #code and size
        if c==0:        #match
            exon_bound.append([ chrom_st,chrom_st + s-1])

            chrom_st += s

        elif c==1:        #insertion to ref
            continue
        elif c==2:        #deletion to ref
            chrom_st += s
        elif c==3:        #gap or intron
            chrom_st += s
        elif c==4:        #soft clipping. We do NOT include soft clip as part of exon
            chrom_st += s
        else:
            continue
    return exon_bound


def ovp_gene(geneIdx,chrom1,exon_blocks1,chrom2,exon_blocks2,strand) :
    
        #(genes1,cds1,utr5_1,utr3_1,intron1,itgUp1k_1,itgDn1k_1) = geneIdx.Gene_annotation(chrom1,exon_blocks1,strand)
        #(genes2,cds2,utr5_2,utr3_2,intron2,itgUp1k_2,itgDn1k_2) = geneIdx.Gene_annotation(chrom2,exon_blocks2,strand)
        
        res1 = geneIdx.Gene_annotation(chrom1,exon_blocks1,strand)
        res2 = geneIdx.Gene_annotation(chrom2,exon_blocks2,strand)
            
        #intersect genes1 and genes2
        gene = ""
        type = ""
        type1 = ""
        type2 = ""
        cds1 = False
        utr5_1 = False
        utr3_1 = False
        intron1 = False
        itgUp1k_1 = False
        itgDn1k_1 = False
        
        
        if len(res1) > 1 and len(res2) > 1 :
            type1 = res1[0]
            type2 = res2[0]
            genes1 = list(set(res1[1:]))
            genes2 = list(set(res2[1:]))
            
            if len(genes1) == len(genes2) == 1:
                if genes1[0] == genes2[0]:
                    gene = genes1[0]
                    if type1 == "cds" or type2 == "cds" :
                        return (gene,True,utr5_1,utr3_1,intron1,itgUp1k_1,itgDn1k_1)
                    if type1 == "utr5" or type2 == "utr5" :
                        return (gene,cds1,True,utr3_1,intron1,itgUp1k_1,itgDn1k_1)
                    if type1 == "utr3" or type2 == "utr3" :
                        return (gene,cds1,utr5_1,True,intron1,itgUp1k_1,itgDn1k_1)
                    if type1 == "intron" or type2 == "intron" :
                        return (gene,cds1,utr5_1,utr3_1,True,itgUp1k_1,itgDn1k_1)
                    if type1 == "Up1k" or type2 == "Up1k" :
                        return (gene,cds1,utr5_1,utr3_1,intron1,True,itgDn1k_1)
                    if type1 == "Dn1k" or type2 == "Dn1k" :
                        return (gene,cds1,utr5_1,utr3_1,intron1,itgUp1k_1,itgDn1k_1)
                else:
                    return (gene,cds1,utr5_1,utr3_1,intron1,itgUp1k_1,itgDn1k_1)
                    
            else:
                return (gene,cds1,utr5_1,utr3_1,intron1,itgUp1k_1,itgDn1k_1)

                   

        if len(res1) == 2  and len(res2) <= 1 :
            gene = res1[1]
            type = res1[0]
        if len(res2) == 2 and len(res1) <= 1 :
            gene = res2[1]
            type = res2[0]

        if type == "cds" :
            cds1 = True
        if type == "utr5" :
            utr5_1 = True
        if type == "utr3" :
            utr3_1 = True

        if type == "intron" :
            intron1 = True
        if type == "Up1k" :
            itgUp1k_1 = True
        if type == "Dn1k" :
            itgDn1k_1 = True

        return (gene,cds1,utr5_1,utr3_1,intron1,itgUp1k_1,itgDn1k_1)

# parsing BAM file
class parseBAM:
    '''This class provides fuctions to parsing SAM or BAM files and quality controls.'''
    
    def __init__(self,inputFile,stranded):
        '''constructor. input could be bam or sam'''
        self.stranded = stranded
        self.max_mapQ = 0
        try:
            self.samfile = pysam.AlignmentFile(inputFile,'rb')
            if len(self.samfile.header) ==0:
                print >>sys.stderr, "BAM/SAM file has no header section. Exit!"
                sys.exit(1)
            self.format = "BAM"
        except:
            self.samfile = pysam.AlignmentFile(inputFile,'r')
            if len(self.samfile.header) ==0:
                print >>sys.stderr, "BAM/SAM file has no header section. Exit!"
                sys.exit(1)
            self.format = "SAM"

#def main(arg
#    rRNAidx = args[1]
#    geneIdx = args[2]
    
    def qc (self,rRNAidx,geneIdx,mapq,outfile_fig,outfile) :

        q_cut = mapq
        
        res = Results()
        '''Calculate mapping statistics'''
        rRNA_read = 0
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
        #multi_reads = []
        multi_read1 = []
        multi_read2 = []
        alignments_per_read = []
        
        clip_prof = Clipping_prof(outfile,outfile_fig,mapq)
        cov_prof = CoverageBody_prof(outfile,outfile_fig,geneIdx)
        rDup_prof = ReadDup_prof(outfile,outfile_fig)
        inDist_prof = InnerDist_prof(outfile,outfile_fig)
        
        
        try:
            while(1):
                #read in one alignment
                aligned_read = self.samfile.next()
                cur_reads = []
                
                #SE reads
                if not aligned_read.is_paired :                    
                    #count unmapped read
                    if aligned_read.is_unmapped :
                        res.unmapped_reads += 1
                        res.total_reads += 1
                        continue
                                    
                    if aligned_read.qname == prev_read_id :
                        alignments_per_read.append((aligned_read,None))
                        continue
                    else :                
                        cur_read = None        
                        if len(alignments_per_read) == 1 : #unique read
                            res.uniq_mapped_reads += 1
                            res.total_reads += 1
                                
                            cur_read = alignments_per_read[0][0]
                            
                            skip_read  = 0
                            if cur_read.is_qcfail or cur_read.mapq < q_cut:    #skip QC fail read
                                clip_prof.set("",cur_read.mapq)
                                res.low_qual += 1
                                skip_read = 1
                            if cur_read.is_duplicate:        #skip duplicate read
                                skip_read = 1
                            if skip_read == 1 :
                                cur_read = None
                            
                        if len(alignments_per_read) > 1: #multi reads
                            res.multi_mapped_reads += 1
                            res.total_reads += 1
                            cur_read = alignments_per_read[0][0]
                            clip_prof.set("",cur_read.mapping_quality)
                            cur_read = None
                        
                        if cur_read is not None :
                            cur_reads.append((cur_read,None))
                            if cur_read.is_reverse:
                                res.reverse_read += 1
                            else:
                                res.forward_read += 1
                                
                        alignments_per_read = []
                        alignments_per_read.append((aligned_read,None))
                        prev_read_id = aligned_read.qname
                        
                #pair end read
                if aligned_read.is_paired :
                    res.is_pairEnd = True 
                    flag_pos = aligned_read.qname.find('/')
                    cur_read_id = aligned_read.qname
                    if flag_pos != -1:
                        cur_read_id = aligned_read.qname[:flag_pos]
                    if cur_read_id == prev_read_id :
                        if aligned_read.is_read1 :
                            multi_read1.append(aligned_read)
                        if aligned_read.is_read2 :
                            multi_read2.append(aligned_read)                        
                    else :
                        cur_read1 = None
                        cur_read2 = None
                        res.total_reads += 1
                        
                        #multi-reads
                        if len(multi_read1) >1 or len(multi_read2) > 1 :
                            res.multi_mapped_reads += 1
                            if len(multi_read1) > 1 :
                                cur_read1 = multi_read1[0]
                                clip_prof.set("",cur_read1.mapping_quality)
                            if len(multi_read2) > 1 :
                                cur_read2 = multi_read2[0]
                                clip_prof.set("",cur_read2.mapping_quality)
                            cur_read1 = None
                            cur_read2 = None
                        else :
                        #uniq read
                            if len(multi_read1) == 1:
                                cur_read1 = multi_read1[0]
                            if len(multi_read2) == 1:
                                cur_read2 = multi_read2[0]

                            mapped = 0
                            if cur_read1 is not None and not cur_read1.is_unmapped :
                                res.mapped_read1 += 1
                            
                                if cur_read1.is_qcfail or cur_read1.mapq < q_cut :
                                    clip_prof.set("",cur_read1.mapping_quality)
                                    cur_read1 = None
                                    res.low_qual_read1 += 1
                                else :
                                    if cur_read1.is_reverse :
                                        res.reverse_read += 1
                                    else:
                                        res.forward_read +=1
                                mapped = 1
                            else :
                                cur_read1 = None
                                if prev_read_id != "" :
                                    res.unmapped_read1 += 1
                        
                            if cur_read2 is not None and not cur_read2.is_unmapped :
                                res.mapped_read2 += 1
                                if cur_read2.is_qcfail or cur_read2.mapq < q_cut :
                                    clip_prof.set("",cur_read2.mapping_quality)
                                    res.low_qual_read2 += 1
                                    cur_read2 = None
                                else :
                                    if cur_read1 is None :
                                        if not cur_read2.is_reverse :
                                            res.reverse_read += 1
                                        else:
                                            res.forward_read +=1
                                mapped = 1
                            else :
                                cur_read2 = None
                                if prev_read_id != "" : #not the first fragment
                                    res.unmapped_read2 += 1
                         
                            if mapped == 1 :
                                res.uniq_mapped_reads += 1
                            else :
                                res.unmapped_reads += 1
                                cur_read1 = None
                                cur_read2 = None
                            
                            if cur_read1 is not None or cur_read2 is not None:
                                cur_reads.append((cur_read1,cur_read2))
                            
                        multi_read1 = []
                        multi_read2 = []
                        prev_read_id = cur_read_id 
                        if aligned_read.is_read1 :
                            multi_read1.append(aligned_read)
                        if aligned_read.is_read2 :
                            multi_read2.append(aligned_read)
                                
                #overlap with genomic features
                for (cur_read1,cur_read2) in cur_reads :
                    strand1 = "."
                    strand2 = "."
                    chrom1 = ""
                    chrom2 = ""
                    exons1 = []
                    exons2 = []
                    RNA_read1 = ""
                    RNA_read2 = ""
                    intron_blocks1 = []
                    intron_blocks2 = []
                    
                    if cur_read1 is not None :
                        RNA_read1 = cur_read1.seq.upper()
                        #mappability computed on read
                        clip_prof.set(cur_read1.cigarstring,cur_read1.mapping_quality)

                        strand1 = "-" if cur_read1.is_reverse else "+"
                        chrom1 = self.samfile.getrname(cur_read1.tid)
                        
                        exons1 = fetch_exon(chrom1, cur_read1.pos, cur_read1.cigar,self.format)
                        #sys.stderr.write(cur_read1.qname+"\t")
                        intron_blocks1 = fetch_intron(cur_read1.pos,cur_read1.cigartuples,self.format)
                        
        
                    if cur_read2 is not None :
                        RNA_read2 = cur_read2.seq.upper()
                        #mappability computed on read
                        #sys.stderr.write(cur_read2.qname+"\n")
                        clip_prof.set(cur_read2.cigarstring,cur_read2.mapping_quality)
                        
                        strand2 = "-" if cur_read2.is_reverse else "+"
                        chrom2 = self.samfile.getrname(cur_read2.tid)
                        
                        exons2 = fetch_exon(chrom2, cur_read2.pos, cur_read2.cigar,self.format)
                        #sys.stderr.write(cur_read2.qname+"\t")
                        intron_blocks2 = fetch_intron(cur_read2.pos,cur_read2.cigartuples,self.format)

                    if len(intron_blocks1) + len(intron_blocks2) == 0:
                        res.noSplice += 1
                    else :
                        res.splice += 1

                    if self.stranded =="reverse" :
                        strand1 = "-" if strand1 == "+" else "+"
                        strand2 = "-" if strand2 == "+" else "+"
                    
                    #paired read
                    if chrom1 != "" and chrom2 != "" :
                        res.paired_reads += 1
                        
                        if  chrom1 != chrom2 :
                            res.paired_diff_chrom += 1
                    
                        if strand1 == "-" and strand2 == "-" :
                            res.mapped_minus_minus  += 1
                        if strand1 == "+" and strand2 == "+" :
                            res.mapped_plus_plus += 1
                        if strand1 == "+" and strand2 == "-" :
                            res.mapped_plus_minus += 1
                        if strand1 == "-" and strand2 == "+" :
                            res.mapped_minus_plus += 1
                                
                    # non-stranded
                    if self.stranded == "no" :
                        strand1 = "."
                        strand2 = "."

                    if cur_read1 is None :
                        strand1 = "+" if strand2 == "-" else "-"

                    #read duplicate rate
                    rDup_prof.count(RNA_read1,RNA_read2,chrom1,chrom2,exons1,exons2,strand1,strand2)
                    #rRNA contamination
                    if rRNAidx is not None :
                        if rRNAidx.is_rRNA(chrom1,exons1,strand1) or rRNAidx.is_rRNA(chrom2,exons2,strand1):
                            rRNA_read += 1
                            continue
                    cur_time = time.time()
                    #print("time before ovp_gene " + str(cur_time))
                    (gene,cds,utr5,utr3,intron,itgUp1k,itgDn1k) = ovp_gene(geneIdx,chrom1,exons1,chrom2,exons2,strand1)
                    #cur_time = time.time()
                    cov_prof.count(gene,exons1,exons2)
                    if cur_read1 is not None and cur_read2 is not None:
                        #estimate insertion size
                        if cur_read1.is_proper_pair :
                            if strand1 == "+" or strand1 == "." :
                                if self.stranded == "reverse" :
                                    inDist_prof.count(geneIdx,gene,cur_read2,chrom2,intron_blocks2,strand1)
                                else:
								    inDist_prof.count(geneIdx,gene,cur_read1,chrom1,intron_blocks1,strand1)
                            else :
                                if self.stranded == "reverse" :
								    inDist_prof.count(geneIdx,gene,cur_read1,chrom1,intron_blocks1,strand1)
                                else :
                                    inDist_prof.count(geneIdx,gene,cur_read2,chrom2,intron_blocks2,strand1)
                    

                    if cds :
                        cds_exon_read += 1
                        continue

                    if utr5 :
                        utr_5_read += 1
                        continue

                    if utr3 :
                        utr_3_read += 1
                        continue

                    if intron :
                        intron_read += 1
                        continue

                    if itgUp1k :
                        intergenic_up1kb_read += 1
                        continue

                    if itgDn1k :
                        intergenic_down1kb_read += 1
                        continue

                    intergenic_read += 1
                    
                                
        except StopIteration:
            pass
        
        res.splice = int(res.splice)
        res.noSplice = int(res.noSplice)
        
        if res.is_pairEnd :
            clip_prof.write(res.total_reads*2-res.unmapped_read1-res.unmapped_read2)
        else :
            clip_prof.write(res.total_reads)
        rDup_prof.write()
        inDist_prof.write()
        (total_exons,zero_exons) = cov_prof.write(res.total_reads)

        res.read_dist_plot_file1 = outfile_fig + ".read_distr.png"
        res.read_dist_plot_file2 = outfile_fig + ".read_distr_pie.png"
        
        try :
            ROUT = open(outfile + '.read_distr.r', 'w')
            print >> ROUT, 'png("' + outfile_fig + '.read_distr.png",width=500,height=500,units="px")\n'
            print >> ROUT, "M=c(" + str(cds_exon_read) + "," + str(utr_5_read) + "," + str(utr_3_read) + "," + str(intron_read) + "," + str(intergenic_up1kb_read) + "," + str(intergenic_down1kb_read) + "," + str(rRNA_read) + "," + str(intergenic_read) + ")\n"
            
            print >> ROUT, 'Mname=c("CDS","5UTR","3UTR","Intron","TSS_Up_1Kb","TES_Down_1Kb","rRNA","Others")\n'
            print >> ROUT, 'val = barplot(M,xlab="",space=1,ylab="Read Counts",col="blue",border="NA")\n'
            print >> ROUT, 'text(x=seq(val[1],val[8],by=2),y=rep(0,8),srt=60,adj=0,offset=2,pos=1,xpd=T,labels=Mname)\n'
            print >> ROUT, "dev.state = dev.off()"
            ROUT.close()

            ROUT = open(outfile + '.read_distr_pie.r', 'w')
            if total_exons != 0 :
                print >> ROUT,'png("' +outfile_fig + '.read_distr_pie.png",width=500,height=500,units="px")\n'
                print >> ROUT, 'pie(c(' + str(total_exons-zero_exons) + ',' + str(zero_exons) + '),labels=c("Covered","Uncovered"),main="Exons",radius=0.6,clockwise=T,col=c("blue","white"))\n'
                print >> ROUT, "dev.state = dev.off()"
            ROUT.close()
        except :
            sys.stderr.write("Error in writing plotting scripts.\n")
            pass

        
        res.insert_plot_file = inDist_prof.plot_file
        res.insert_file = inDist_prof.out_file2
        res.clipping_plot_file = clip_prof.plot_file
        res.mapq_plot_file = clip_prof.plot_file2
        res.mapq_file = clip_prof._out_file3
        res.read_dup_plot_file = rDup_prof.plot_file
        res.readLen_plot_file = rDup_prof.plot_file2
        res.read_cov_plot_file = cov_prof.plot_file
        res.geneCount_file = cov_prof.outfile3
        res.trans_cov_plot_file = cov_prof.plot_file2
        res.seqDeDup_percent = rDup_prof.seqDeDup_percent
        res.posDeDup_percent = rDup_prof.posDeDup_percent
        
        return res

#if __name__ == '__main__':
#    main(args)



