//
//  InnerDist_prof.cpp
//  BAMQC_c++
//
//  Created by Ying Jin on 10/28/15.
//  Copyright (c) 2015 Ying Jin. All rights reserved.
//

#include "InnerDist_prof.h"
#include "Mappability.h"
#include <map>
#include <iostream>
#include <fstream>
//#include <stdio>
#include <stdlib.h>
#include <math.h>
#include <htslib/sam.h>

InnerDist_prof::InnerDist_prof(std::string outfile_data,std::string outfile_fig,int sample_size,int l_bound,int u_bound,int s)
{
    
    InnDist_data_file = outfile_data + ".inner_distance_freq.txt";
    InnDist_script_file= outfile_data + ".inner_distance_plot.r";
    InnDist_fig_file = outfile_fig+".inner_distance_plot.png";
    lower_bound = l_bound;
    upper_bound = u_bound;
    
    //samplesize = sample_size;
    step = s;
    
    //window_left_bound = range(lower_bound,up_bound,step);
    
}
InnerDist_prof::InnerDist_prof(int sample_size,int l_bound,int u_bound,int s)
{
    lower_bound = l_bound;
    upper_bound = u_bound;
    
    //samplesize = sample_size;
    step = s;
}
InnerDist_prof::~InnerDist_prof(){
    
}
void InnerDist_prof::write()
{
    //estimate the inner distance of mRNA pair end fragment. fragment size = insert_size + 2 x read_length'''
    if (pair_num == 0){
        return ;
    }
    try{

        std::string pos_str = "";
        std::string cnt_str = "";
        int size = ceil((upper_bound-lower_bound)/step)-1;
        int to_plot[size];
        
        for(int i=0;i<size;i++)
        {
            to_plot[i] = 0;
        }
        for (auto& kv: counts){
            if (kv.first >= lower_bound && kv.first <= upper_bound) {
                int pos = ceil((kv.first - lower_bound - step/2)/step);
                if (pos < 0) {
                    pos = 0;
                }
                to_plot[pos] += kv.second;
            }
        }
        std::ofstream FQ ;
        FQ.open(InnDist_data_file,std::ofstream::out);

        for(int i=0; i < size; i++) {
            int pos = (int)(i*step + lower_bound + step/2);
            int pos1 = i*step + lower_bound;
            
            pos_str +=','+ std::to_string(pos);
            cnt_str += ',' + std::to_string(to_plot[i]);
            
            FQ << pos1 << '\t' << std::to_string(pos1 + step-1) << '\t' << to_plot[i] << '\n' ;
        }
		FQ.close();
        if(pos_str != ""){
        std::ofstream RS ;
        RS.open(InnDist_script_file,std::ofstream::out);
        //plot_file = outfile_fig+".inner_distance_plot.png"
        RS << "png(\"" << InnDist_fig_file << "\",width=500,height=500,units=\"px\")\n";

        RS << "fragsize=rep(c(" << pos_str.substr(1) << "),times=c(" << cnt_str.substr(1) << + "))\n";
        RS << "frag_sd = round(sd(fragsize),digits=0)\n";
        RS << "frag_mean = round(mean(fragsize),digits=0)\n";
        RS << "frag_median = round(median(fragsize),digits=0)\n";
        
        RS << "hist(fragsize,probability=T,breaks="<< size + 1 <<",xlim=c("<< lower_bound << ',' <<upper_bound << "),xlab=\"Inner distance (bp)\",main=paste(c(\"Median=\",frag_median,\";\",\"Mean=\",frag_mean,\";\",\"SD=\",frag_sd),collapse=\"\"),border=\"blue\")\n";
        
        RS << "lines(density(fragsize,bw=" << 2 * step <<"),col='red')\n";
        RS << "abline(v=frag_median,lty=2)\n";
        RS << "dev.state = dev.off()\n";
        
        RS.close();
		}
    }catch(std::ofstream::failure e){
        std::cout << "Error in writing inner distance profile." << std::endl;
        return;
    }
}
        
void InnerDist_prof::add(InnerDist_prof * inDist)
{
    //samplesize += inDist->samplesize;
    pair_num += inDist->pair_num;

    for(auto& kv: inDist->counts)
    {
        std::map<int, int>::iterator it;
        it = counts.find(kv.first);
        if (it != counts.end()) {
            counts[kv.first] += kv.second;
        }
        else {
            counts[kv.first] = kv.second;
        }
    }
}

void InnerDist_prof::count(GeneFeatures * geneIdx,int type,bam1_t * aligned_read,std::string chrom,std::vector<std::pair<int,int> > intron_blocks,std::string strand)
{
    //if (pair_num >= samplesize)
    //    return ;
    int splice_intron_size=0;
    //int read1_len = aligned_read->core.l_qseq; //infer_query_length()
    uint32_t *cigar = bam_get_cigar(aligned_read);
    int read1_len = bam_cigar2mapped_read_len(aligned_read->core.n_cigar,cigar);
    
    int read1_start = aligned_read->core.pos;
    int read2_start = aligned_read->core.mpos;
    int read1_end = 0;
    
    for(auto& intron : intron_blocks){
        splice_intron_size += intron.second - intron.first;
    }
    
    read1_end = read1_start + read1_len + splice_intron_size;
    int inner_distance = read2_start - read1_end +1;
  //  if(inner_distance > -130 && inner_distance < -120)
//{  
  //          char * name = bam_get_qname(aligned_read);
//	printf("read name %s \n", name);
    
//}
    if (inner_distance >= lower_bound && inner_distance <= upper_bound)
    {
        pair_num ++;
        std::vector<std::pair<int,int> > exons ;
        if (type == CDS || type== UTR5 || type== UTR3 ) {
            exons = geneIdx->get_exons(chrom,read1_end,read2_start,strand);
        }
        
        if (exons.size() > 0 )
        {
            int size = 0;
            for (auto& p : exons){
                size += p.second - p.first;
            }
            if (size <= inner_distance && size > 1){
                counts[size] += 1;
            }
            else{
                counts[inner_distance] += 1;
            }
        }
        else{
            counts[inner_distance] += 1;
        }
    }
            
}
