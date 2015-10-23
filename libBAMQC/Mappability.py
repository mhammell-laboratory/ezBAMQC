#!/usr/bin/env python
'''mappability.'''

#import built-in modules
import os,sys,string, re
import collections
import math


        
class Clipping_prof:
    def __init__(self,outfile_data,outfile_fig,qcut):
        self._out_file1 = outfile_data + ".clipping_profile.xls"
        self._out_file2 = outfile_data + ".clipping_profile.r"
        self.plot_file = outfile_fig+ ".clipping_profile.png"
        self.plot_file2 = outfile_fig+ ".mapq_profile.png"
        self._out_file3 = outfile_data + ".mapq_profile.xls"
        self._out_file4 = outfile_data + ".mapq_profile.r"
        self.soft_p = re.compile(r'(.*?)(\d+)S')
        self.read_part = re.compile(r'(\d+)[MIS=X]')
        #total_read =0
        #self.skip_part_of_read =0
        self.soft_clip_profile=collections.defaultdict(int)
        self.mapqlist = collections.defaultdict(int)
        self.q_cutoff = qcut
        #self.read_pos=[]
        #self.clip_count=[]
    
    def write(self,total_read):
        OUT=open(self._out_file1,'w')
        read_pos = []
        clip_count = []
        
        print >>OUT, "Position\tRead_Total\tRead_clipped"
        if len(self.soft_clip_profile) > 0 :
            max_i = max(self.soft_clip_profile.keys())
            for j in range(max_i+1) :
                if j not in self.soft_clip_profile :
                    self.soft_clip_profile[j] = 0
            self.soft_clip_profile[0] = 1
            for i in self.soft_clip_profile:
                read_pos.append(str(i))
                clip_count.append(str(self.soft_clip_profile[i]))
            
                print >>OUT, str(i) + '\t' + str(total_read) + '\t' + str(self.soft_clip_profile[i])

        OUT.close()
        
        OUT=open(self._out_file3,'w')
        mapq_val = []
        mapq_count = []
        
        print >>OUT, "MAPQ\tRead_Total\tRead_with_mapq"
                
        max_mapq = 0
        for q in self.mapqlist:
            mapq_val.append(str(q))
            mapq_count.append(str(self.mapqlist[q]))
            if q > max_mapq :
                max_mapq = q
            print >>OUT, str(q) + '\t' + str(total_read) + '\t' + str(self.mapqlist[q])
        
        OUT.close()
            
        if len(self.soft_clip_profile) > 0 :
            try :
                ROUT=open(self._out_file2,'w')
                
                print >>ROUT, "png(\"%s\",width=500,height=500,units='px')\n" % (self.plot_file)
                print >>ROUT, "read_pos=c(" + ','.join(read_pos) + ')\n'
                print >>ROUT, 'read_pos=read_pos[2:length(read_pos)]\n'
                print >>ROUT, "count=c(" + ','.join(clip_count) + ')\n'
                print >>ROUT, "count=count[2:length(count)]\n"
                
                print >>ROUT, 'plot(read_pos,1-(count/'+str(total_read)+'),pch=20,xlab="Position of reads",ylab="Mappability",col="blue")\n'
                print >>ROUT, "dev.state=dev.off()\n"
                ROUT.close()
            except :
                sys.stderr.write("Error in writing clipping plot script.\n")
                pass
                    
        if len(self.mapqlist) > 0 :
            try :
                ROUT=open(self._out_file4,'w')

                print >>ROUT, "png(\"%s\",width=500,height=500,units='px')\n" % (self.plot_file2)
                print >>ROUT, "mapq_val=c(" + ','.join(mapq_val) + ')\n'
                print >>ROUT, "mapq_count=c(" + ','.join(mapq_count) + ')\n'
                
                print >>ROUT, 'xname=c("<3","<10","<20","<30","30-' +str(max_mapq)+'")'
                print >>ROUT, 'freq = rep(0,5)'
                
                print >>ROUT, 'freq[1] = sum(mapq_count[which(mapq_val<3)])/' + str(total_read) + '*100'
                print >>ROUT, 'freq[2] = sum(mapq_count[which(mapq_val<10)])/' + str(total_read) + '*100'
                print >>ROUT, 'freq[3] = sum(mapq_count[which(mapq_val<20)])/' + str(total_read) + '*100'
                print >>ROUT, 'freq[4] = sum(mapq_count[which(mapq_val<30)])/' + str(total_read) + '*100'
                print >>ROUT, 'freq[5] = 100'
                
                print >>ROUT, 'barplot(freq,beside=T,xlab="Mapping Quality",border="NA",space=1.5,main="Mapping Quality",ylim=c(0,100),ylab="Cumulative proportion (%)",col="blue",names.arg=xname)'
                #print >>ROUT, 'abline(v=' + str(self.q_cutoff) + ', lty=2)\n'
                print >>ROUT, "dev.state=dev.off()\n"
                ROUT.close()
            except :
                sys.stderr.write("Error in writing mapq plot script.\n")
                pass

    
    def set(self,cigar_str,mapq):
        
        self.mapqlist[mapq] += 1
        m = self.soft_p.findall(cigar_str)
        
        skip_part_of_read =1
        if len(m) > 0 :
            for j in m: 
                skip_part_of_read += sum([int(i) for i in self.read_part.findall(j[0])])
                for n in range(skip_part_of_read,(skip_part_of_read + int(j[1])+1)):
                    self.soft_clip_profile[n]+=1

                skip_part_of_read += int(j[1])



