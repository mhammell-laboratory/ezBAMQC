//
//  InnerDist_prof.h
//  BAMQC_c++
//
//  Created by Ying Jin on 10/28/15.
//  Copyright (c) 2015 Ying Jin. All rights reserved.
//

#ifndef __BAMQC_c____InnerDist_prof__
#define __BAMQC_c____InnerDist_prof__

#include <stdio.h>
#include <vector>
#include <utility>

#include "GeneFeatures.h"
#include "Constants.h"
#include "htslib/sam.h"

class InnerDist_prof{
public:
    
    std::string InnDist_data_file;
    std::string InnDist_fig_file;
    std::string InnDist_script_file;
    int samplesize = 0;
    int step= 0;
    int lower_bound = 0;
    int upper_bound = 0;
    int pair_num = 0;
    
    std::map<int,int> counts ;
    
    InnerDist_prof(int sample_size=SAMPLESIZE,int low_bound= LOW_BOUND,int up_bound= UPPER_BOUND,int step= STEP);
    InnerDist_prof(std::string data_dir,std::string fig_dir,int sample_size=SAMPLESIZE,int low_bound= LOW_BOUND,int up_bound= UPPER_BOUND,int step= STEP);
    ~InnerDist_prof();
    void write();
    void count(GeneFeatures * geneIdx, int type,bam1_t * aligned_read,std::string chrom,std::vector<std::pair<int,int> > intron_blocks,std::string strand);
    void add(InnerDist_prof * inDist_prof);
};

#endif /* defined(__BAMQC_c____InnerDist_prof__) */
