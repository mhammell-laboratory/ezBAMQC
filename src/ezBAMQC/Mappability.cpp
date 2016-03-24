//
//  Mappablity.cpp
//  BAMQC_c++
//
//  Created by Ying Jin on 10/28/15.
//  Copyright (c) 2015 Ying Jin. All rights reserved.
//

#include "Mappability.h"
#include "htslib/sam.h"
//#include <stdio>
#include <iostream>
#include <fstream>

#include <stdlib.h>
#include <vector>


int bam_cigar2mapped_read_len(int n_cigar, const uint32_t *cigar)
{
    int k, l;
	    for (k = l = 0; k < n_cigar; ++k)
		{
                  if (bam_cigar_op(cigar[k])==BAM_CMATCH) 
				  {          l += bam_cigar_oplen(cigar[k]);}
	    }
		   return l;
}

void bam_cigar2Clip(int n_cigar, const uint32_t *cigar,std::vector<std::pair<int, int> > * clip_pos_len)
{
    int k, pos;

    //std::cout << "number of cigar " << n_cigar <<std::endl;  

    for (k = pos = 0; k < n_cigar; ++k)
        if (bam_cigar_op(cigar[k]) == BAM_CSOFT_CLIP )
        {
            clip_pos_len->push_back(std::pair<int,int> (pos,bam_cigar_oplen(cigar[k])));
            pos += bam_cigar_oplen(cigar[k]);
        }
        else {
		if (bam_cigar_op(cigar[k]) == BAM_CMATCH || bam_cigar_op(cigar[k]) == BAM_CINS || bam_cigar_op(cigar[k]) == BAM_CEQUAL || bam_cigar_op(cigar[k]) == BAM_CDIFF )
            { pos += bam_cigar_oplen(cigar[k]);}
        }

    return ;
}


Clipping_prof::Clipping_prof(std::string outfile_data,std::string outfile_fig,int qcut,int sample_size)
{
    
    
    clip_data_file = outfile_data + ".clipping_profile.xls";
    clip_script_file = outfile_data + ".clipping_profile.r";
    clip_fig_file = outfile_fig+ ".clipping_profile.png";
    
    mapq_fig_file = outfile_fig+ ".mapq_profile.png";
    mapq_data_file = outfile_data + ".mapq_profile.xls";
    mapq_script_file = outfile_data + ".mapq_profile.r";
    
    readlen_data_file = outfile_data + ".readlen_profile.xls";
    readlen_fig_file = outfile_fig+ ".readlen_profile.png";
    readlen_script_file = outfile_data + ".ReadLen_plot.r";
    
    samplesize = sample_size;
    
    q_cutoff = qcut;
    read_len =0;
}
Clipping_prof::~Clipping_prof(){
    
}
Clipping_prof::Clipping_prof(int mapq,int sample_size){
    q_cutoff = mapq;
    samplesize = sample_size;
    read_len =0;

}
int Clipping_prof::get_max_read_len(){
    return this->read_len;
}
void Clipping_prof::write(int total_read)
{
    std::string read_pos ="";
    std::string clip_count = "";
    std::string mapq_val = "";
    std::string mapq_count = "";
    std::string readlen_val = "";
    std::string readlen_count = "";
    
    int max_mapq = 0;
    
    if (read_len > MAX_READ_LEN) {
        std::cout << "read length greater than 10000." << std::endl;
        return;
    }
    try{
        std::ofstream OUT ;
        OUT.open (readlen_data_file, std::ofstream::out);
        OUT << "Position\tRead_Total\tRead_Len_mapped\n";
        
        //soft_clip_profile[0] = 1;
        for(auto& kv : readLen_list ){
            readlen_val += ',' + std::to_string(kv.first);
            readlen_count += ',' + std::to_string(kv.second);
            
            OUT << kv.first << '\t' << total_read << '\t' << kv.second << '\n';
        }
        
        OUT.close();
        
        
    }catch(std::ofstream::failure e ){
        std::cout << "Error in writing clipping profile." << std::endl;
        return;
    }
    
    try{
        std::ofstream OUT ;
        OUT.open (clip_data_file, std::ofstream::out);
        OUT << "Position\tRead_Total\tRead_clipped\n";
        
        //soft_clip_profile[0] = 1;
        for(int i=0; i< read_len; i++ ){
            read_pos += ',' + std::to_string(i);
            clip_count += ',' + std::to_string(soft_clip_profile[i]);
        
            OUT << i << '\t' << total_read << '\t' << soft_clip_profile[i] << '\n';
        }

        OUT.close();
        
        
    }catch(std::ofstream::failure e ){
        std::cout << "Error in writing clipping profile." << std::endl;
        return;
    }
    try{
	    if(read_pos != ""){
        std::ofstream ROUT ;
        ROUT.open(clip_script_file,std::ofstream::out);
        
        ROUT << "png(\"" << clip_fig_file << "\",width=500,height=500,units=\"px\")\n";
        ROUT << "read_pos=c(" << read_pos.substr(1) << ")\n" ;
        //ROUT <<'read_pos=read_pos[2:length(read_pos)]\n';
        ROUT << "count=c(" << clip_count.substr(1) << ")\n";
        //ROUT << "count=count[2:length(count)]\n";
                
        ROUT << "plot(read_pos,1-(count/" << total_read <<"),pch=20,xlab=\"Position of reads\",ylab=\"Mappability\",col=\"blue\")\n";
        
        ROUT << "dev.state=dev.off()\n";
        
        ROUT.close();
		}
    }catch(std::ofstream::failure e ){
        std::cout << "Error in writing clipping script file." << std::endl;
        return;
    }
    try{
        std::ofstream OUT;
        OUT.open(mapq_data_file,std::ofstream::out);
        
        OUT << "MAPQ\tRead_Total\tRead_with_mapq\n";
        
        for(auto& kv :mapqlist ){
            mapq_val += ',' + std::to_string(kv.first);
            
            mapq_count += ',' + std::to_string(kv.second);
            
            if (kv.first > max_mapq ){
                max_mapq = kv.first;
            }
            
            OUT << std::to_string(kv.first) +  '\t' + std::to_string(total_read) + '\t' + std::to_string(kv.second) +'\n';
        }
        OUT.close();
    }catch(std::ofstream::failure e ){
        std::cout << "Error in writing mapping quality profile." << std::endl;
        return;
    }

    if (mapqlist.size() > 0 ){
        try {
		if(mapq_val != ""){
            std::ofstream ROUT;
            ROUT.open(mapq_script_file,std::ofstream::out);

            ROUT << "png(\"" << mapq_fig_file << "\",width=500,height=500,units=\"px\")\n";
            ROUT << "mapq_val=c(" <<  mapq_val.substr(1) << ")\n";
            ROUT<< "mapq_count=c(" << mapq_count.substr(1) << ")\n";
                
            ROUT << "xname=c(\"<3\",\"<10\",\"<20\",\"<30\",\"30-"  << std::to_string(max_mapq) << "\")\n";
            ROUT <<"freq = rep(0,5)\n";
                
            ROUT << "freq[1] = sum(mapq_count[which(mapq_val<3)])/" << std::to_string(total_read) << "*100\n";
            ROUT <<"freq[2] = sum(mapq_count[which(mapq_val<10)])/" << std::to_string(total_read) << "*100\n";
            ROUT << "freq[3] = sum(mapq_count[which(mapq_val<20)])/" << std::to_string(total_read) << "*100\n";
            ROUT << "freq[4] = sum(mapq_count[which(mapq_val<30)])/" << std::to_string(total_read) << "*100\n";
            ROUT<< "freq[5] = 100\n";
                
            ROUT << "barplot(freq,beside=T,xlab=\"Mapping Quality\",border=\"NA\",space=1.5,main=\"Mapping Quality\",ylim=c(0,100),ylab=\"Cumulative proportion (%)\",col=\"blue\",names.arg=xname)\n";
            
            ROUT << "dev.state=dev.off()\n";
            ROUT.close();
            }
        }catch(std::ofstream::failure e ){
            std::cout << "Error in writing mapping quality script file." << std::endl;
            return;
        }
    }
	//std::cout << readLen_list.size() << std::endl;
    if (readLen_list.size() > 0) {
        try {
            std::ofstream ROUT;
            ROUT.open(readlen_script_file,std::ofstream::out);
            
            ROUT << "png(\"" << readlen_fig_file << "\",width=500,height=500,units=\"px\")\n";
            ROUT << "readlen_val=c(" <<  readlen_val.substr(1) << ")\n";
            ROUT<< "readlen_count=c(" << readlen_count.substr(1) << ")\n";
            
            ROUT << "plot(readlen_val,(readlen_count/" << total_read <<"),pch=20,xlab=\"Mapped Read Length\",ylab=\"Proportion\",col=\"blue\")\n";
            
            ROUT << "dev.state=dev.off()\n";
            ROUT.close();
            
        }catch(std::ofstream::failure e ){
            std::cout << "Error in writing mapping quality script file." << std::endl;
            return;
        }
    }
}

void Clipping_prof::add(Clipping_prof * clip_prof)
{
    for (auto& j : clip_prof->mapqlist){
        std::map<int,int>::iterator it;
        it = mapqlist.find(j.first);
        
        if (it != mapqlist.end()){
            mapqlist[j.first] += j.second;
        }
        else {
            mapqlist[j.first] = j.second;
        }
    }
    
    for (int i=0;i< MAX_READ_LEN; i++){
        soft_clip_profile[i] += clip_prof->soft_clip_profile[i];
    }
    if (clip_prof->get_max_read_len() > read_len) {
        read_len = clip_prof->get_max_read_len();
    }
    for (auto& kv : clip_prof->readLen_list){
        std::map<int,int>::iterator it;
        it = readLen_list.find(kv.first);
        
        if (it != readLen_list.end()){
            readLen_list[kv.first] += kv.second;
        }
        else {
            readLen_list[kv.first] = kv.second;
        }
    }
}

void Clipping_prof::set_qual(int mapq)
{
    std::map<int,int>::iterator it;
    it = mapqlist.find(mapq);
    if (it != mapqlist.end()) {
        mapqlist[mapq] += 1;
    }
    else{
        mapqlist[mapq] = 1;
    }
}
void Clipping_prof::set(int len,uint32_t n_cigar, uint32_t * cigar,int mapq)
{
    int r = bam_cigar2mapped_read_len(n_cigar,cigar);
	//std::cout << "read len " << r << std::endl;
    //std::cout << "number of cigar " << n_cigar <<std::endl;  
    if(read_len < len)
    {
        read_len = len;
    }
    std::map<int,int>::iterator it;
    it = mapqlist.find(mapq);
    if (it != mapqlist.end()) {
        mapqlist[mapq] += 1;
    }
    else{
        mapqlist[mapq] = 1;
    }
	//std::cout << "after mapq " << mapq << std::endl;
    it = readLen_list.find(r);
    if(it!=readLen_list.end()){
        readLen_list[r] += 1;
    }
    else{
        readLen_list[r] = 1;
    }
    
	//std::cout << "after read len " << r << std::endl;
    std::vector<std::pair<int, int> > soft_clip_pos_len ;
    
    bam_cigar2Clip(n_cigar,cigar,&soft_clip_pos_len);
    
    
	//std::cout << "soft clip len " << soft_clip_pos_len.size() << std::endl;
    if(soft_clip_pos_len.size() > 0) {
        for (auto p : soft_clip_pos_len){
            int pos = p.first;
            int len = p.second;
            
	//std::cout << pos  << "\t" << len << std::endl;
            for (int j = 0; j< len;j++){
                soft_clip_profile[pos+j] +=1;
            }
            
        }
    }

}

