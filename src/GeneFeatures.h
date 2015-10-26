//
//  GeneFeatures.h
//  BAMQC-0.5
//
//  Created by Ying Jin on 9/15/15.
//  Copyright (c) 2015 Ying Jin. All rights reserved.
//

#ifndef __BAMQC_0_5__GeneFeatures__
#define __BAMQC_0_5__GeneFeatures__

#include <stdio.h>
#include <string.h>
#include <vector>
#include <map>

#include "IntervalTree.h"

#define BIN_SIZE 1000000
//typedef struct ITV{
 //   long start;
  //  long end;
//} ITV;

typedef struct {
    std::string chrom;
    long start;
    long end;
} chr_ITV;

extern "C" bool itv_comp(Interval, Interval);



typedef std::map<std::string, std::vector<std::string> > gene_exon_Dict ;
typedef gene_exon_Dict::iterator gene_exon_Dict_It;

typedef std::map<std::string,std::map<int, IntervalTree *> > chrom_itvTree_Dict ;
typedef chrom_itvTree_Dict::iterator chrom_itvTree_Dict_itr;
//typedef Dict::const_iterator gene_exon_Dict_It;


class Gene {
public:
    std::vector<std::pair<long,long> > exons;
    std::vector<std::pair<long,long> > cds;
    std::string id;
    std::string strand;
    std::vector<std::pair<long,long> > utr5;
    std::vector<std::pair<long,long> > utr3 ;
    std::vector<std::pair<long,long> > intron ;
    std::vector<std::pair<long,long> > itg1k ;
    Gene(std::string gid, std::string ss);
    //Gene& operator= (const Gene& gg);
    //Gene(const std::string gid, const std::string ss);
    ~Gene();
    void add_cds(long st, long end);
    void add_exons(long st,long end);
    void get_others();
};

class GeneFeatures{

public:
    std::vector<std::string> features;
    
    gene_exon_Dict plus_gene_exons ;
    gene_exon_Dict minus_gene_exons ;
    
    //chrom_itvTree_Dict plus_features;
    //chrom_itvTree_Dict minus_features;
    
    chrom_itvTree_Dict cds_exon_idx_plus ;
    chrom_itvTree_Dict cds_exon_idx_minus ;
    
    GeneFeatures(std::string GTFfilename,std::string id_attribute);
    ~GeneFeatures();
    //void set(std::string GTFfilename,std::string id_attribute);
    void read_features(std::string gff_filename, std::string id_attribute) ;
    std::vector<std::pair<long,long> > find_exons(std::string chrom,long st,long end,std::string strand,std::string gene) ;
    std::vector<std::string> getFeatures() ;
    std::vector<std::string> get_exon_pair(std::string g);
    //std::vector<std::string> ovp_others(std::vector<chr_ITV> itv_list,std::string strand) ;
    //std::vector<std::string> Gene_annotation(std::vector<chr_ITV> itv_list, std::string strand);
    std::vector<std::string> Gene_annotation(std::string, std::vector<std::pair<long,long> > itv_list, std::string strand);
private:
    void build_tree(std::map<std::string, std::map<std::string,Gene> > temp_plus, std::map<std::string, std::map<std::string,Gene> > temp_minus);

};

/*extern "C" {
    void quick_sort(std::vector<Interval> &intervals, int first, int last);
    int pivot(std::vector<Interval> &intervals, int first, int last);
    long get_first(const std::pair<long, long>& p);
    
    long get_last(const std::pair<long, long>& p);
};*/

#endif /* defined(__BAMQC_0_5__GeneFeatures__) */
