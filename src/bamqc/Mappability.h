//
//  Mappablity.h
//  BAMQC_c++
//
//  Created by Ying Jin on 10/28/15.
//  Copyright (c) 2015 Ying Jin. All rights reserved.
//

#ifndef __BAMQC_c____Mappablity__
#define __BAMQC_c____Mappablity__

#include <stdio.h>
#include <string>
#include "Constants.h"
#include <utility>
#include <map>

int bam_cigar2mapped_read_len(int n_cigar, const uint32_t *cigar);
class Clipping_prof{
public:
    int q_cutoff;
    int samplesize;
    
    int soft_clip_profile[MAX_READ_LEN]={0};
    std::string clip_data_file;
    std::string clip_fig_file;
    std::string clip_script_file;
    
    std::string readlen_data_file;
    std::string readlen_script_file;
    std::string readlen_fig_file;
    
    std::string mapq_fig_file;
    std::string mapq_data_file;
    std::string mapq_script_file;
    
    std::map<int, int> mapqlist;
    std::map<int, int> readLen_list;
    
    Clipping_prof(int qcut, int sample_size = SAMPLESIZE);
    Clipping_prof(std::string data_dir,std::string fig_dir,int qcut, int sample_size = SAMPLESIZE);
    
    ~Clipping_prof();
    void set(int len,uint32_t n_cigar, uint32_t * cigar,int mapq);
    void set_qual(int mapq);
    void write(int total_read);
    void add(Clipping_prof * clip_prof);

    int get_max_read_len();
private:
    int read_len;
};

#endif /* defined(__BAMQC_c____Mappablity__) */
