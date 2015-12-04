//
//  Results.cpp
//  BAMQC_c++
//
//  Created by Ying Jin on 10/28/15.
//  Copyright (c) 2015 Ying Jin. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>

#include "Results.h"

Results::Results(){
    
}

Results::~Results(){
    
}

void Results::write(std::string fname)
{
    try{
        std::ofstream f;
        f.open(fname,std::ofstream::out);
        
        f << "filename" << "\t" << fname << "\n";
        f << "is_pairEnd" << "\t" << this->is_pairEnd << "\n";
        f << "clipping_plot_file" << "\t" << this->clipping_plot_file << "\n";
        f << "mapq_plot_file"<< "\t" <<this->mapq_plot_file << "\n";
        f << "mapq_file"<< "\t" <<this->mapq_file << "\n";
        f << "read_cov_plot_file"<< "\t" <<this->read_cov_plot_file << "\n";
        f << "trans_cov_plot_file"<< "\t" <<this->trans_cov_plot_file<< "\n";
        f << "insert_plot_file" << "\t" << this->insert_plot_file<< "\n";
        f << "insert_file" << "\t" << this->insert_file<< "\n";
        f <<  "read_dist_plot_file1" << "\t" << this->read_dist_plot_file1<< "\n";
        f <<  "read_dist_plot_file2" << "\t" <<this->read_dist_plot_file2<< "\n";
        f <<  "readLen_plot_file" << "\t" << this->readLen_plot_file<< "\n";
        f <<  "geneCount_file" << "\t" << this->geneCount_file<< "\n";
        f <<  "seqDeDup_percent" << "\t" <<std::to_string(this->seqDeDup_percent)<< "\n";
        f <<  "posDeDup_percent" << "\t" << std::to_string(this->posDeDup_percent)<< "\n";
        f <<  "no_clipping" << "\t" << this->no_clipping << "\n";
        f <<  "no_rRNA" << "\t" << this->no_rRNA << "\n";
        f <<  "total_reads" << "\t" <<std::to_string(this->total_reads) << "\n";
        f <<  "uniq_mapped_reads" << "\t" <<std::to_string(this->uniq_mapped_reads)<< "\n";
        f <<  "multi_mapped_reads" << "\t" <<std::to_string(this->multi_mapped_reads)<< "\n";
        f <<  "unmapped_reads" << "\t" <<std::to_string(this->unmapped_reads)<< "\n";
        f << "low_qual" << "\t" <<std::to_string(this->low_qual)<< "\n";
        f << "low_qual_read1" << "\t" <<std::to_string(this->low_qual_read1)<< "\n";
        f << "low_qual_read2" << "\t" <<std::to_string(this->low_qual_read2)<< "\n";
        f << "pcr_dup" << "\t" <<std::to_string(this->pcr_dup)<< "\n";
        f << "rRNA_read" << "\t" <<std::to_string(this->rRNA_read)<< "\n";
        
        f << "cds_read" << "\t" <<std::to_string(this->cds_exon_read)<< "\n";
        f << "utr5_read" << "\t" <<std::to_string(this->utr_5_read)<< "\n";
        f << "utr3_read" << "\t" <<std::to_string(this->utr_3_read)<< "\n";
        f << "intron_read" << "\t" <<std::to_string(this->intron_read)<< "\n";
        f << "itgup1k_read" << "\t" <<std::to_string(this->intergenic_up1kb_read)<< "\n";
        f << "itgdn1k_read" << "\t" <<std::to_string(this->intergenic_down1kb_read)<< "\n";
        f << "itg_read" << "\t" <<std::to_string(this->intergenic_read)<< "\n";
        
        
        f << "unmapped_read1" << "\t" <<std::to_string(this->unmapped_read1)<< "\n";
        f << "unmapped_read2" << "\t" <<std::to_string(this->unmapped_read2)<< "\n";
        f << "mapped_read1"<< "\t" <<std::to_string(this->mapped_read1)<< "\n";
        f << "mapped_read2"<< "\t" <<std::to_string(this->mapped_read2) << "\n";
        f << "forward_read"<< "\t" <<std::to_string(this->forward_read)<< "\n";
        f << "reverse_read" << "\t" << std::to_string(this->reverse_read)<< "\n";
        f << "paired_reads"<< "\t" << std::to_string(this->paired_reads) + "\n";
        f << "mapped_plus_minus"<< "\t" <<std::to_string(this->mapped_plus_minus) << "\n";
        f << "mapped_plus_plus"<< "\t" <<std::to_string(this->mapped_plus_plus)<< "\n";
        f << "mapped_minus_plus" << "\t" <<std::to_string(this->mapped_minus_plus)<< "\n";
        f << "mapped_minus_minus" << "\t" <<std::to_string(this->mapped_minus_minus)<< "\n";
        f << "ins_read" << "\t" <<std::to_string(this->ins_read)<< "\n";
        f << "del_read"<< "\t" <<std::to_string(this->del_read)<< "\n";
        f << "noSplice"<< "\t" <<std::to_string(this->noSplice)<< "\n";
        f << "splice"<< "\t" <<std::to_string(this->splice)<< "\n";
        f << "paired_diff_chrom" << "\t" <<std::to_string(this->paired_diff_chrom)<< "\n";
        f.close();

    }catch(std::ofstream::failure e ){
        std::cout << "Error in writing result ." << std::endl;
        return;
    }

}

void Results::add(Results *res){
    
    no_clipping = (no_clipping || res->no_clipping);
    no_rRNA = (no_rRNA && res->no_rRNA);
    
    rRNA_read += res->rRNA_read;
    total_reads += res->total_reads;
    uniq_mapped_reads +=res->uniq_mapped_reads;
    multi_mapped_reads +=res->multi_mapped_reads;
    unmapped_reads += res->unmapped_reads;
    low_qual += res->low_qual;
    low_qual_read1 += res->low_qual_read1;
    low_qual_read2 += res->low_qual_read2;
    pcr_dup += res->pcr_dup;
    
    unmapped_read1 += res->unmapped_read1;
    unmapped_read2 += res->unmapped_read2;
    mapped_read1 += res->mapped_read1;
    mapped_read2 += res->mapped_read2;
    forward_read += res->forward_read;
    reverse_read += res->reverse_read;
    paired_reads += res->paired_reads;
    
    mapped_plus_minus += res->mapped_plus_minus;
    mapped_plus_plus += res->mapped_plus_plus;
    mapped_minus_plus += res->mapped_minus_plus;
    mapped_minus_minus += res->mapped_minus_minus;
    
    ins_read += res->ins_read;
    del_read += res->del_read;
    
    noSplice += res->noSplice;
    splice += res->splice;
    paired_diff_chrom += res->paired_diff_chrom;

    cds_exon_read += res->cds_exon_read;
    utr_5_read += res->utr_5_read;
    utr_3_read += res->utr_3_read;
    intron_read += res->intron_read;
    intergenic_up1kb_read += res->intergenic_up1kb_read;
    intergenic_down1kb_read += res->intergenic_down1kb_read;
    intergenic_read += res->intergenic_read;
}
