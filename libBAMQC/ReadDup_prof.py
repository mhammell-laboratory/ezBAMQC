#!/usr/bin/env python
'''PCR duplicates.'''

#import built-in modules
import os,sys
import re
import string
import collections


#changes to the paths
class ReadDup_prof:
    def __init__(self,outfile_data,outfile_fig):
        self.outfile1 = outfile_data+".seq.DupRate.xls"
        self.outfile2 = outfile_data + ".pos.DupRate.xls"
        self.outfile3 = outfile_data + ".DupRate_plot.r"
        self.outfile4 = outfile_data + ".ReadLen_plot.r"
        self.outfile5 = outfile_data + ".ReadLen.xls"
        self.total = 0
        self.total_readLen = 0
        self.seqDup=collections.defaultdict(int)
        self.posDup=collections.defaultdict(int)
        self.readLenList = collections.defaultdict(int)
        self.plot_file = outfile_fig+".DupRate_plot.png"
        self.plot_file2 = outfile_fig + ".ReadLength_plot.png"
        self.seqDeDup_percent = 0
        self.posDeDup_percent = 0
        
        
    def write(self,up_bound=500):
        '''Calculate reads's duplicate rates'''
        SEQ=open(self.outfile1,'w')
        POS=open(self.outfile2,'w')
        RS=open(self.outfile3,'w')
        RL = open(self.outfile5,'w')
        RL_plot = open(self.outfile4,'w')

        seqDup_count=collections.defaultdict(int)
        posDup_count=collections.defaultdict(int)


        print >>RL, "ReadLength\tFrequency"
        for key in self.readLenList.keys() :
            print >>RL, str(key) + "\t" + str(round(1.0*self.readLenList[key]/(self.total_readLen),4))
        RL.close()

        print >>RL_plot, "png(\'%s\',width=500,height=500,units='px')\n" % (self.plot_file2)
        
        print >>RL_plot, "seq_occ=c(" + ','.join([str(i) for i in sorted(self.readLenList.keys()) ]) + ')'
        #print >>RL_plot, "seq_freq=c(" + ','.join([str(round(1.0*self.readLenList[i]/(self.total_readLen),4)) for i in sorted(self.readLenList.keys()) ]) + ')'
        print >>RL_plot, "seq_freq=c(" + ','.join([str(round(100.0*self.readLenList[i]/(self.total_readLen),4)) for i in sorted(self.readLenList.keys()) ]) + ')'
        print >>RL_plot, "seq_perc=rep(0,length(seq_freq))"
        print >>RL_plot, "for(i in 1:length(seq_occ)) {seq_perc[i] = sum(seq_freq[which(seq_occ>=seq_occ[i])])}"

        print >>RL_plot, "plot(x=seq_occ,y=seq_perc, ylab='Proportion of reads >= the read length (%)',xlab='Read Length',pch=20,cex=0.8,col='blue',type='b')"
        print >>RL_plot, 'dev.state = dev.off()'
        RL_plot.close()

        print >>SEQ, "Occurrence\tUniqReadNumber"
        
        seqDeDup_cnt = 0
        for i in self.seqDup.values():    #key is occurence, value is uniq reads number (based on seq)
            seqDup_count[i] +=1
        
        for k in sorted(seqDup_count.iterkeys()):
            print >>SEQ, str(k) +'\t'+ str(seqDup_count[k])
            seqDeDup_cnt += seqDup_count[k]
        SEQ.close()
        
        #sys.stderr.write("total : "+str(self.total)+"\n")
        #sys.stderr.write(str(seqDeDup_cnt)+"\n")

        self.seqDeDup_percent = round((1.0*seqDeDup_cnt/self.total),4)*100
        
        
        posDeDup_cnt = 0
        print >>POS, "Occurrence\tUniqReadNumber"
        for i in self.posDup.values():    #key is occurence, value is uniq reads number (based on coord)
            posDup_count[i] +=1
        
        for k in sorted(posDup_count.iterkeys()):
            print >>POS, str(k) +'\t'+ str(posDup_count[k])
            posDeDup_cnt += posDup_count[k]
        POS.close()
        self.posDeDup_percent = round((1.0*posDeDup_cnt/self.total),4)*100

        #sys.stderr.write(str(self.seqDeDup_percent)+"\n")
        #     sys.stderr.write(str(self.posDeDup_percent)+"\n")

        title ="Reads remaining if deduplicated by sequence " + str(self.seqDeDup_percent)+"%" + "\n" + "Reads remaining if deduplicated by position " + str(self.posDeDup_percent) + "%"
        print >>RS, "png(\'%s\',width=500,height=500,units='px')\n" % (self.plot_file)
        print >>RS, 'main = "' + title + '"\n'
        print >>RS, "seq_occ=c(" + ','.join([str(i) for i in sorted(seqDup_count.iterkeys()) ]) + ')'
        print >>RS, "seq_uniqRead=c(" + ','.join([str(seqDup_count[i]) for i in sorted(seqDup_count.iterkeys()) ]) + ')'
        print >>RS, "pos_occ=c(" + ','.join([str(i) for i in sorted(posDup_count.iterkeys()) ]) + ')'
        print >>RS, "pos_uniqRead=c(" + ','.join([str(posDup_count[i]) for i in sorted(posDup_count.iterkeys()) ]) + ')'
        print >>RS," total = sum(pos_uniqRead*pos_occ)"
        print >>RS,"pos_occ_perc= log2(pos_uniqRead/total*100)"
        print >>RS,"seq_occ_perc= log2(seq_uniqRead/total*100)"
        print >>RS,"ymax = "+str(up_bound)
        print >>RS,"if(max(c(seq_occ,pos_occ)) + 100 < ymax) {ymax = max(c(seq_occ,pos_occ)) + 100}"
        print >>RS, "plot(pos_occ,pos_occ_perc,pch=20,col='blue',yaxt='n',xlim=c(1,"+ str(up_bound)+"),ylab='Percentage (%)',xlab='#Duplicates',main=main)"
        print >>RS, 'ym=ceiling(max(pos_occ_perc))'
        print >>RS, 'ymin = floor(min(pos_occ_perc))'
        print >>RS, 'lab = round(2^seq(ymin,ym,by=2),2)'
        print >>RS, 'lab[which(lab==0)]=""'
        print >>RS, 'axis(2,at=seq(ymin,ym,by=2),labels=lab)'
        print >>RS, 'points(seq_occ,seq_occ_perc,pch=20,col="red")'

        print >>RS, "legend(%d,ym,bty='n',legend=c('Sequence-base','Mapping-base'),col=c('red','blue'),pch=c(20,20))\n" % max(up_bound-200,1)

        print >>RS, 'dev.state = dev.off()'
        
        RS.close()
    
    def count(self,RNA_read1,RNA_read2,chrom1,chrom2,exon_blocks1,exon_blocks2,strand1,strand2):
                exon_boundary=""
                self.total += 1
                seq = RNA_read1 + "_" + RNA_read2
                read1_len = len(RNA_read1)
                read2_len = len(RNA_read2)
                if read1_len > 0 :
                    self.readLenList[read1_len] += 1
                    self.total_readLen += 1
                if read2_len > 0 :
                    self.readLenList[read2_len] += 1
                    self.total_readLen += 1
                try :
                    self.seqDup[seq] +=1                    #key is read sequence
                except :
                    self.seqDup[seq] = 1

                exons1 = sorted(exon_blocks1,key=lambda tup:tup[0])
                exons2 = sorted(exon_blocks2,key=lambda tup:tup[0])
                key = ""

                if len(exon_blocks1) == 0 and len(exon_blocks2) > 0 :
                    key = chrom2 + ":" + str(exons2[0][0]) + ":" + strand2

                if len(exon_blocks2)  == 0 and len(exon_blocks1) > 0 :
                    key = chrom1 + ":" + str(exons1[0][0]) + ":" + strand1

                if len(exon_blocks1) > 0 and len(exon_blocks2) > 0  :
                    chr1 = chrom1 #exons1[0][0]
                    chr2 = chrom2 #exons2[0][0]
                    
                    key = chr1 + ":" +chr2 + ":" + strand1 + ":" + strand2 + ":"
                    if strand1 == "+" or strand1 == "." :
                        key += str(exons1[0][0]) + ":" + str(exons2[-1][1])
                    if strand1 == "-" :
                        key += str(exons2[0][0]) + ":" + str(exons1[-1][1])
                if key is not "" :
                    try:
                        self.posDup[key] +=1
                    except :
                        self.posDup[key] = 1





