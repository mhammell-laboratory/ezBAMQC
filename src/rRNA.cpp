//
//  GeneFeatures.cpp
//  BAMQC-0.5
//
//  Created by Ying Jin on 9/15/15.
//  Copyright (c) 2015 Ying Jin. All rights reserved.
//

#include "rRNA.h"

#include <cmath>
#include <fstream>
#include <sstream>
//#include <regex>
#include "stdlib.h"
#include <algorithm>
#include <iostream>



//long get_first(const std::pair<long, long>& p) { return p.first; }

//long get_last(const std::pair<long, long>& p) { return p.second; }

//bool itv_comp(Interval first, Interval second){
//    return first.start < second.start ;
//}

//int pivot(std::vector<Interval> &intervals, int first, int last)
/*{
    int  p = first;
    long pivotElement = intervals[first].start;
    
    
    for(int i = first+1 ; i <= last ; i++)
    {
        // If you want to sort the list in the other order, change "<=" to ">"
        if(intervals[i].start <= pivotElement)
        {
            std::swap(intervals[i],intervals[p]);
            p++;
            
        }
    }
    
    return p;
}*/

/*void quick_sort(std::vector<Interval> &intervals, int first, int last){
    
    int pivotElement;
    
    if(first < last)
    {
        pivotElement = pivot(intervals, first, last);
        quick_sort(intervals, first, pivotElement-1);
        quick_sort(intervals, pivotElement+1, last);
    }
    
}*/


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
    long start;
    long end;
    
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
                    (*tmp_ptr).push_back(Interval(feature,start,end,"rRNA"));
                }
                else{
                    std::vector<Interval> tmp ;
                    tmp.push_back(Interval(feature,start,end,"rRNA"));

                    temp_plus.insert(std::pair<std::string,std::vector<Interval> >(chrom,tmp));
                }
            }
                    
            if (strand == "-" ) {
                
                tmp_itr = temp_minus.find(chrom);
                if (tmp_itr != temp_minus.end()) {
                    std::vector<Interval> *tmp_ptr = &(tmp_itr->second);
                    (*tmp_ptr).push_back(Interval(feature,start,end,"rRNA"));
                }
                else{
                    std::vector<Interval> tmp ;
                    tmp.push_back(Interval(feature,start,end,"rRNA"));
                    
                    temp_minus.insert(std::pair<std::string,std::vector<Interval> >(chrom,tmp));
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

    
        //trees for CDS
        //sys.stderr.write("start to build interval tree\n")
    //std::cout << "temp_plus size " << temp_plus.size() << std::endl;
    //std::cout << "temp_plus chrom1 size" << temp_plus["chr1"].size() << std::endl;
    //build_tree(temp_plus,temp_minus);
   // std::cout << "after built tree" << std::endl;
    /*for (tmp_itr = temp_plus.begin(); tmp_itr != temp_plus.end(); tmp_itr++) {
        gene_id_map = tmp_itr->second;
        for (id_itr = gene_id_map.begin(); id_itr != gene_id_map.end(); id_itr ++) {
            delete id_itr->second;
        }
        delete gene_id_map;
    }
    for (tmp_itr = temp_minus.begin(); tmp_itr != temp_minus.end(); tmp_itr++) {
        gene_id_map = tmp_itr->second;
        for (id_itr = gene_id_map.begin(); id_itr != gene_id_map.end(); id_itr ++) {
            delete id_itr->second;
        }
        delete gene_id_map;
    }*/
    
}

//find exons of given gene that overlap with the given intervals
//return list of tuples
bool rRNA::is_rRNA(std::string chrom, std::vector<std::pair<long,long> > itv_list,std::string strand)
{
    std::vector<std::string> rRNAs;
    std::vector<Interval> fs ;
   // std::vector<Interval> tmp ;
    std::map<std::string, IntervalTree*>::iterator chrom_it;
    //std::map<int,IntervalTree *>::iterator bin_iter;
    //int i,bin_id_s, bin_id_e;
    int i;
    
    if (strand == "+" || strand == ".") {
        chrom_it = featureIdxs_plus.find(chrom);
        if (chrom_it !=  featureIdxs_plus.end()) {
            
            for (i=0; i < itv_list.size(); i ++) {
                    //std::cout << "start to search tree" << std::endl;
                    fs = (*featureIdxs_plus[chrom]).find(itv_list[i].first,itv_list[i].second);
                                        //std::cout << fs.size() << std::endl;
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

                                        //std::cout << fs.size() << std::endl;
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
