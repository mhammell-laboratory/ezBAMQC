//
//  Coverage_prof.cpp
//  BAMQC_c++
//
//  Created by Ying Jin on 10/28/15.
//  Copyright (c) 2015 Ying Jin. All rights reserved.
//

#include "Coverage_prof.h"

#include <string>
//#include <stdio>
#include <iostream>
#include <fstream>
#include <stdlib.h>



Coverage_prof::Coverage_prof(std::string outfile_data,std::string outfile_fig, GeneFeatures * geneIdx)
{
    frag_num = 0;
    cov_script_file = outfile_data + ".geneBodyCoverage_plot.r";
    cov_data_file = outfile_data + ".geneBodyCoverage.txt";
    cov_fig_file = outfile_fig + ".geneBodyCoverage.png";
    
    transcov_fig_file = outfile_fig + ".TransCoverage.png";
    transcov_script_file = outfile_data + ".TransCoverage.r";
    transcov_data_file = outfile_data + ".geneAbundance.txt";
    
    int num_of_genes = geneIdx->get_numofgenes();
    int num_of_exons = geneIdx->total_exon;
    //geneCounts_list[num_of_genes];
    for (int i=0; i< num_of_genes; i++) {
        geneCounts_list.push_back(0);
		std::vector<int> idx_cnt ;
		for(int j=0;j<=100; j++)
		{
		   idx_cnt.push_back(0);
		}
		gene_percentile_base.push_back(idx_cnt);
    }
	//std::cout << num_of_genes << std::endl;
    for (int i=0; i< num_of_exons; i++) {
        mapped_exon.push_back(0);
    }
    gene_Idx = geneIdx;
    total_exons = num_of_exons;
}

Coverage_prof::Coverage_prof(GeneFeatures * geneIdx){
    gene_Idx = geneIdx;
    frag_num = 0;
    total_exons = geneIdx->total_exon;
    
    int num_of_genes = geneIdx->get_numofgenes();
    int num_of_exons = geneIdx->total_exon;
    //geneCounts_list[num_of_genes];
    for (int i=0; i< num_of_genes; i++) {
        geneCounts_list.push_back(0);
		std::vector<int> idx_cnt ;
		for(int j=0;j<=100; j++)
		{
		   idx_cnt.push_back(0);
		}
		gene_percentile_base.push_back(idx_cnt);
    }
	//std::cout << num_of_genes << std::endl;
    for (int i=0; i< num_of_exons; i++) {
        mapped_exon.push_back(0);
    }

}
    
int Coverage_prof::write(int totalReads){
    int coverage[101] = {0};
    //total_exons = 0;
    int zero_exons = 0;
    //std::cout << "start to write" << std::endl;

    for(size_t i=0;i< gene_percentile_base.size();i++)
    {
        std::vector<int> percentile_list = gene_percentile_base[i];
        for(size_t j=0;j< percentile_list.size();j++){
            coverage[j] += percentile_list[j];
        
        }
        
    }
    std::string geneCnt_str = "";
    for (size_t i=0;i<geneCounts_list.size(); i++) {
        geneCnt_str += ',' + std::to_string(geneCounts_list[i]);
    }
    std::string x_coord="";
    std::string y_coord="";
    try {
        std::ofstream OUT1, OUT2, OUT3, OUT4;
        OUT2.open(cov_data_file,std::ofstream::out);
        OUT3.open(transcov_data_file,std::ofstream::out);
		OUT3 << "gene" << "\t" << "Counts\n";
        if (geneCnt_str !=""){
		OUT4.open(transcov_script_file,std::ofstream::out);
        
        OUT4 << "png(\'" << transcov_fig_file << "\',width=500,height=500,units='px')\n";
        OUT4 << "a=c("<< geneCnt_str.substr(1) << ")\n";
        OUT4 << "Fn = ecdf(a)\n";
        OUT4 << "max_x = round(log(max(knots(Fn)),2),0)\n";
        OUT4 << "xx = c(0,2^seq(0,max_x,by=2))\n";
        OUT4 << "y=Fn(xx)\n";
        OUT4 << "xlog = log(xx[2:length(xx)],base=2)\n";
        OUT4 << "plot(x=c(-1,xlog),y=y,xaxt = 'n',type=\"b\",col=\"blue\",pch=20,xlab=\"Number of Reads\",ylab=\"Cumulative proportion of Genes\")\n";
        OUT4 << "axis(1,at = c(-1,seq(0,max_x,by=2)),labels=c(0,2^seq(0,max_x,by=2)))\n";
        OUT4 << "dev.state = dev.off()";
        OUT4.close();
        }    
        OUT2 << "Total reads: " << std::to_string(totalReads) << "\n";
        OUT2 << "Fragment number: " << frag_num << "\n";
        OUT2 << "percentile\tcount\n";
            
        for(int i=0; i< 101;i++){
            x_coord += ',' + std::to_string(i);
            y_coord += ',' + std::to_string(coverage[i]);
            
            OUT2 << std::to_string(i) << '\t' << std::to_string(coverage[i]) << "\n";
        }
       if(x_coord != ""){ 
        OUT1.open(cov_script_file,std::ofstream::out);
        OUT1 << "png(\'" << cov_fig_file << "\',width=500,height=500,units='px')\n";
        OUT1 << "x=c(" << x_coord.substr(1) << ")\n";
        OUT1 << "y=c(" << y_coord.substr(1) << ")\n";
        
        OUT1 << "smoothsp = smooth.spline(x,y,spar=0.35)\n";
        OUT1 << "plot(smoothsp,type=\"l\",col=\"blue\",xlab=\"Percentile of Gene Body (5\'->3\')\",ylab=\"Number of read\",xlim=c(0,100))\n";
        OUT1 << "dev.state = dev.off()";
            
        OUT1.close();
		}
        OUT2.close();
        
        for(size_t i=0;i<geneCounts_list.size();i++)
        {
            
            OUT3 << gene_Idx->get_name(i) << "\t" << std::to_string(geneCounts_list[i]) << "\n";
        }
        OUT3.close();
    
    }catch(std::ofstream::failure e ){
        std::cout << "Error in writing clipping profile." << std::endl;
        return -1;
    }
    zero_exons = total_exons;
    for (size_t i=0; i< mapped_exon.size(); i++) {
        if (mapped_exon[i] > 0) {
            zero_exons --;
        }
    }
    return zero_exons;
}

void Coverage_prof::add(Coverage_prof * cov_prof)
{
    frag_num += cov_prof->frag_num;
    
    for(size_t i=0;i< geneCounts_list.size();i++){
        geneCounts_list[i] += cov_prof->geneCounts_list[i];
    }
    
    for(size_t i=0;i<gene_percentile_base.size();i++){
        std::vector<int> * g0_perc_cnt = &(gene_percentile_base[i]);
        std::vector<int> g1_perc_cnt = cov_prof->gene_percentile_base[i];
        for (size_t j=0; j<= 100; j++ ) {
            g0_perc_cnt->at(j) = g1_perc_cnt[j] + g0_perc_cnt->at(j);
        }
    }
    for (int i=0; i< total_exons; i++) {
        mapped_exon[i] += cov_prof->mapped_exon[i];
    }
}

void Coverage_prof::count(int gene,std::vector<std::pair<int,int> > exon_blocks1,std::vector<std::pair<int,int> > exon_blocks2,std::vector<int> exons)
{
    for (size_t i=0; i < exons.size(); i++) {
        mapped_exon[exons[i]] =1;
    }
    if (gene != -1){
        frag_num += exon_blocks1.size() + exon_blocks2.size();
        geneCounts_list[gene] += 1;
        
        //__per_base_count(gene,exon_blocks,exon_blocks2);
        int gene_start = gene_Idx->get_start(gene);
        int gene_stop  = gene_Idx->get_stop(gene);
        for (auto& kv : exon_blocks1) {
            if (kv.second >= gene_start && kv.first <= gene_stop)
            {
                for (int j=kv.first; j<=kv.second; j++) {
                    int idx = gene_Idx->exist_in_percentile_list(gene,j);
                    if (idx != -1) {
                        gene_percentile_base[gene][idx] +=1;
                    }
                }
            }
        }
        
        for (auto& kv : exon_blocks2) {
            if (kv.second >= gene_start && kv.first <= gene_stop)
            {
                for (int j=kv.first; j<=kv.second; j++) {
                    int idx = gene_Idx->exist_in_percentile_list(gene,j);
                    if (idx != -1) {
                        gene_percentile_base[gene][idx] +=1;
                    }
                }
            }
        }
        
    }
}

