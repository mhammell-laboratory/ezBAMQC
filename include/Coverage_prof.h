//
//  Coverage_prof.h
//  BAMQC_c++
//
//  Created by Ying Jin on 11/18/15.
//  Copyright (c) 2015 Ying Jin. All rights reserved.
//

#ifndef __BAMQC_c____Coverage_prof__
#define __BAMQC_c____Coverage_prof__

#include <stdio.h>
#include <vector>
#include <utility>

#include "GeneFeatures.h"

class Coverage_prof{
public:
    int frag_num ;
    int total_exons;
    GeneFeatures * gene_Idx;
    
    std::string cov_script_file ;
    std::string cov_data_file ;
    std::string cov_fig_file ;
    
    std::vector<int> geneCounts_list;
    std::vector<int> mapped_exon;
    std::vector<std::vector<int> > gene_percentile_base;
    
    std::string transcov_fig_file ;
    std::string transcov_script_file ;
    std::string transcov_data_file ;
    
    Coverage_prof(GeneFeatures * geneIdx);
    Coverage_prof(std::string outfile_data,std::string outfile_fig, GeneFeatures * geneIdx);
    int write(int totalReads);
    void add(Coverage_prof * cov_prof);
    void count(int gene,std::vector<std::pair<int,int> > exon_blocks1,std::vector<std::pair<int,int> > exon_blocks2,std::vector<int> exons);
};

#endif /* defined(__BAMQC_c____Coverage_prof__) */
