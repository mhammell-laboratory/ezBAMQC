//
//  rRNA.h
//  BAMQC-0.5
//
//  Created by Ying Jin on 9/15/15.
//  Copyright (c) 2015 Ying Jin. All rights reserved.
//

#ifndef __BAMQC_0_5__rRNA__
#define __BAMQC_0_5__rRNA__

#include <stdio.h>
#include <string.h>
#include <vector>
#include <map>
#include <utility>

#include "IntervalTree.h"
#include "GeneFeatures.h"


class rRNA{

public:
    
    std::map<std::string, IntervalTree *>  featureIdxs_plus ;
    std::map<std::string, IntervalTree *>  featureIdxs_minus ;
    
    rRNA(std::string rRNAfilename);
    ~rRNA();
    
    bool is_rRNA(std::string chrom,std::vector<std::pair<int,int> > itv_list,std::string strand);
    
private:
    void read_features(std::string rRNA_filename) ;
    //void build_tree(std::map<std::string, std::map<std::string,Gene> > temp_plus, std::map<std::string, std::map<std::string,Gene> > temp_minus);

};


#endif 
