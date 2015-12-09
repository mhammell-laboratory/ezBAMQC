//
//  GeneFeatures.cpp
//  BAMQC-0.5
//
//  Created by Ying Jin on 9/15/15.
//  Copyright (c) 2015 Ying Jin. All rights reserved.
//

#include "rRNA.h"
#include "Constants.h"

#include <cmath>
#include <fstream>
#include <sstream>
//#include <regex>
#include "stdlib.h"
#include <algorithm>
#include <iostream>


rRNA::rRNA(std::string rRNAfilename){

        read_features(rRNAfilename);

}

rRNA::~rRNA(){
    std::map<std::string,IntervalTree *>::iterator it;
    
    for (it=featureIdxs_plus.begin(); it != featureIdxs_plus.end(); it++) {
        IntervalTree *tmp = it->second;
        delete tmp;
        
    }
    for (it=featureIdxs_minus.begin(); it != featureIdxs_minus.end(); it++) {
        IntervalTree *tmp = it->second;
        delete tmp;
    }
    
}

//Reading & processing annotation files
void rRNA::read_features(std::string rRNA_filename)
{

    //dict of dicts since the builtin type doesn't support it for some reason
    std::map<std::string, std::vector<Interval> > temp_plus ;
    std::map<std::string, std::vector<Interval> > temp_minus ;
    std::map<std::string, std::vector<Interval> >::iterator tmp_itr;

    int i =  0;
    //int counts = 0 ;
    int line_no = 0;
    int start;
    int end;
    int idx = 0;
    
    std::ifstream input; //(gff_filename);
    
    try{
        input.open (rRNA_filename, std::ifstream::in);

        while(! input.eof()){
    
            std::string line,chrom,feature,start_ss,end_ss,score,strand;
            std::stringstream ss;

            if (! std::getline(input,line)){ break; }
            line_no ++;
        
            if (line == "\n" || !line.compare(0,1,"#")) { continue; }
       
            ss << line;
            std::getline(ss,chrom,'\t');
            std::getline(ss,start_ss,'\t');
            std::getline(ss,end_ss,'\t');
            std::getline(ss,feature,'\t');
            std::getline(ss,score,'\t');
            std::getline(ss,strand,'\t');
       
            try{
                start = std::stol(start_ss);
                end = std::stol(end_ss);
            }
            catch (const std::invalid_argument& ia) {
                std::cerr << "Invalid argument: " << ia.what() << '\n';
                std::exit(1);
           
            }
        
            if (strand == "+" ){
                tmp_itr = temp_plus.find(chrom);
                if (tmp_itr != temp_plus.end()) {
                    std::vector<Interval> *tmp_ptr = &(tmp_itr->second);
                    //(*tmp_ptr).push_back(Interval(feature,-1,start,end,"rRNA"));
                    (*tmp_ptr).push_back(Interval(idx,-1,start,end,RRNA));
                    idx++;
                }
                else{
                    std::vector<Interval> tmp ;
                    tmp.push_back(Interval(idx,-1,start,end,RRNA));

                    temp_plus.insert(std::pair<std::string,std::vector<Interval> >(chrom,tmp));
                    idx++;
                }
            }
                    
            if (strand == "-" ) {
                
                tmp_itr = temp_minus.find(chrom);
                if (tmp_itr != temp_minus.end()) {
                    std::vector<Interval> *tmp_ptr = &(tmp_itr->second);
                    (*tmp_ptr).push_back(Interval(idx,-1,start,end,RRNA));
                    idx++;
                }
                else{
                    std::vector<Interval> tmp ;
                    tmp.push_back(Interval(idx,-1,start,end,RRNA));
                    
                    temp_minus.insert(std::pair<std::string,std::vector<Interval> >(chrom,tmp));
                    idx++;
                }
            }
    
            i += 1 ;
            if (i % 100000 == 0 )
            {
            //sys.stderr.write("%d GTF lines processed.\n" % i);
                std::cout << i << " rRNA lines processed." << std::endl;
            }
        
        }
    
        input.close();
    }
    catch(std::ifstream::failure e){
        std::cout << "error in read file " << rRNA_filename << std::endl;
    }
        //build interval trees
    
    
    for (tmp_itr = temp_plus.begin(); tmp_itr != temp_plus.end(); tmp_itr++) {
        std::string chr = tmp_itr->first;
        std::vector<Interval> itemlist = tmp_itr->second;
        
        std::sort(itemlist.begin(),itemlist.end(),itv_comp) ;
        
        featureIdxs_plus[chr] = new IntervalTree(itemlist);
    }
    
    for (tmp_itr = temp_minus.begin(); tmp_itr != temp_minus.end(); tmp_itr++) {
        std::string chr = tmp_itr->first;
        
        std::vector<Interval> itemlist = tmp_itr->second;
        
        std::sort(itemlist.begin(),itemlist.end(),itv_comp) ;
        
        featureIdxs_minus[chr] = new IntervalTree(itemlist);
    }

    
}

//find exons of given gene that overlap with the given intervals
//return list of tuples
bool rRNA::is_rRNA(std::string chrom, std::vector<std::pair<int,int> > itv_list,std::string strand)
{
    std::vector<std::string> rRNAs;
    std::vector<Interval> fs ;
    std::map<std::string, IntervalTree*>::iterator chrom_it;
    size_t i;
    if (strand == "+" || strand == ".") {
        chrom_it = featureIdxs_plus.find(chrom);
        if (chrom_it !=  featureIdxs_plus.end()) {
            
            for (i=0; i < itv_list.size(); i ++) {
                    //std::cout << "start to search tree" << std::endl;
                    fs = (*featureIdxs_plus[chrom]).find(itv_list[i].first,itv_list[i].second);
                //}
            }
        }
    }

        if (strand == "-" or strand == ".") {
            chrom_it = featureIdxs_minus.find(chrom);
            if (chrom_it !=  featureIdxs_minus.end()) {
                
                for (i=0; i < itv_list.size(); i ++) {
                    std::vector<Interval>  tmp = (*featureIdxs_minus[chrom]).find(itv_list[i].first,itv_list[i].second);
                    
                    fs.insert(fs.end(),tmp.begin(),tmp.end());
                }
                
           }
        }
    
    //std::cout << fs.size() << std::endl;
    if (fs.size() > 0 ) {
        return true;
    }
    else{
        return false;
    }
    
}
/*
int main() {
    std::string filename = "./test.bed";
    std::string id_attr = "gene_id";
   
    std::vector<chr_ITV> itv_list ;
    chr_ITV exp ;
    exp.chrom = "chr1";
    exp.start = 11870;
    exp.end = 73220;
    itv_list.push_back(exp);
    
    std::cout << "start to build tree " << std::endl;
    rRNA gIdx (filename);
    std::cout << "after  build tree " << std::endl;

   // for (int i=0; i < itv_list.size(); i++) {
    //    std::cout << itv_list[i].start << std::endl;
    //}
    bool res = gIdx.is_rRNA(itv_list,".");
    //for (int i=0;i<res.size(); i++) {
        std::cout << res << std::endl;
    //}
    
}*/
