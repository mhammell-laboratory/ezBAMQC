//
//  Results.h
//  BAMQC_c++
//
//  Created by Ying Jin on 10/28/15.
//  Copyright (c) 2015 Ying Jin. All rights reserved.
//

#ifndef BAMQC_c___Results_h
#define BAMQC_c___Results_h

#include <string>


class Results {

public:
    Results();
    ~Results();
    void add(Results *res);
    void write(std::string fname);
    
    std::string filename = "";
    
    std::string clipping_plot_file = "";
    std::string mapq_plot_file = "";
    std::string mapq_file = "";
    std::string read_cov_plot_file = "";;
    std::string trans_cov_plot_file = "";
    std::string insert_plot_file = "";
    std::string insert_file = "";
    std::string read_dist_plot_file1 = "";
    std::string read_dist_plot_file2 = "";
    std::string read_dup_plot_file = "";
    std::string readLen_plot_file = "";
    std::string geneCount_file = "";

    double seqDeDup_percent = 0;
    double posDeDup_percent = 0;
    bool is_pairEnd = false;
    
    bool no_clipping = false;
    bool no_rRNA = false;


    int total_reads = 0;
    int uniq_mapped_reads = 0;
    int multi_mapped_reads = 0;
    int unmapped_reads = 0;
    int low_qual = 0;
    int low_qual_read1 = 0;
    int low_qual_read2 = 0;
    int pcr_dup = 0;
    int rRNA_read = 0;
    
    int cds_exon_read =0;
    int utr_5_read =0;
    int utr_3_read =0;
    int intron_read =0;
    int intergenic_up1kb_read =0;
    int intergenic_down1kb_read = 0;
    int intergenic_read = 0;

    int unmapped_read1 = 0;
    int unmapped_read2 = 0;
    int mapped_read1 = 0;
    int mapped_read2 = 0;
    int forward_read = 0;
    int reverse_read = 0;
    int paired_reads = 0;

    int mapped_plus_minus = 0;
    int mapped_plus_plus = 0;
    int mapped_minus_plus = 0;
    int mapped_minus_minus = 0;

    int ins_read = 0;
    int del_read = 0;

    int noSplice = 0;
    int splice = 0;
    int paired_diff_chrom = 0;
    
};


#endif
