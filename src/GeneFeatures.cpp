//
//  GeneFeatures.cpp
//  BAMQC-0.5
//
//  Created by Ying Jin on 9/15/15.
//  Copyright (c) 2015 Ying Jin. All rights reserved.
//

#include "GeneFeatures.h"
#include "IntervalTree.h"
#include <cmath>
#include <fstream>
#include <sstream>
//#include <regex>
#include "stdlib.h"
#include <algorithm>
#include <iostream>
//#include <boost/tokenizer.hpp>

long get_first(const std::pair<long, long>& p) { return p.first; }

long get_last(const std::pair<long, long>& p) { return p.second; }

bool itv_comp(Interval first, Interval second){
    return first.start < second.start ;
}

int pivot(std::vector<Interval> &intervals, int first, int last)
{
    int  p = first;
    long pivotElement = intervals[first].start;
    
    
    for(int i = first+1 ; i <= last ; i++)
    {
        /* If you want to sort the list in the other order, change "<=" to ">" */
        if(intervals[i].start <= pivotElement)
        {
            std::swap(intervals[i],intervals[p]);
            p++;
            
        }
    }
    
    return p;
}

void quick_sort(std::vector<Interval> &intervals, int first, int last){
    
    int pivotElement;
    
    if(first < last)
    {
        pivotElement = pivot(intervals, first, last);
        quick_sort(intervals, first, pivotElement-1);
        quick_sort(intervals, pivotElement+1, last);
    }
    
}



//Gene::Gene(const std::string gid, const std::string ss){
//Gene::Gene(){
//    id = "";
//    strand = "";
//}
Gene::Gene(std::string gid,  std::string ss){
    id = gid;
    strand = ss ;
}
/*Gene& Gene::operator= (const Gene& gg) {
    id = gg.id;
    strand = gg.strand;
    exons = gg.exons;
    cds = gg.cds;
    utr5 = gg.utr5;
    utr3 = gg.utr3;
    intron = gg.intron ;
    itg1k = gg.itg1k;
    return *this;
}*/
Gene::~Gene(){}

void Gene::add_cds(long st, long end){
    std::pair<long,long> cds_interval (st,end);
    cds.push_back(cds_interval);
}

void Gene::add_exons(long st, long end){
    std::pair<long,long> exon_interval (st,end);
    exons.push_back(exon_interval);
}

void Gene::get_others(){
    std::vector<std::pair<long,long> > left_cds;
    std::vector<std::pair<long,long> > left_exons;
    int idx[exons.size()];
    int i,j;
    long itgUp1k_st, itgUp1k_end,itgDn1k_st,itgDn1k_end ;

    sort(exons.begin(),exons.end());
    sort(cds.begin(),cds.end());
    
/*    exon_starts.reserve(std::distance(exons.begin(), exons.end()));
    exon_ends.reserve(std::distance(exons.begin(), exons.end()));
    exon_starts = std::transform(exons.begin(),
                                 exons.end(),
                                 std::back_inserter(exon_starts), get_first);
    
    exon_ends = std::transform(exons.begin(),
                                exons.end(),
                                std::back_inserter(exon_ends), get_last);
*/
    for (i=0;i < cds.size(); i++) {
        for (j=0; j < exons.size(); j++) {
            if (cds[i].first == exons[j].first && cds[i].second == exons[j].second) {
                idx[j] = 1;
                break;
                
            }
        }
        if (j == exons.size()) {
            left_cds.push_back(cds[i]);
        }
    }

    for (i=0; i<exons.size(); i++) {
        if (idx[i] != 1) {
            left_exons.push_back(exons[i]);
        }
    }
    

    for (i = 1; i < exons.size(); i++) {
        intron.push_back(std::make_pair(exons[i-1].second +1, exons[i].first -1));
    }

    if(strand == "+") {
        itgUp1k_st = std::max(long(0),exons[0].first-1000);
        itgUp1k_end = std::max(long(0),exons[0].first-1 );
        itgDn1k_st = exons[exons.size()-1].second + 1;
        itgDn1k_end = exons[exons.size()-1].second + 1000 ;
        itg1k.push_back(std::make_pair(itgUp1k_st,itgUp1k_end)) ;
        itg1k.push_back(std::make_pair(itgDn1k_st,itgDn1k_end)) ;
    

        if (left_cds.size() == 0 ) { utr3 = exons; };
        
        if (left_cds.size() == 1) {
            for( i=0;i < left_exons.size(); i++) {
                
                if (left_exons[i].second < left_cds[0].first) { utr5.push_back(left_exons[i]); }
                if (left_exons[i].first < left_cds[0].first && left_exons[i].second > left_cds[0].first ) { utr5.push_back(std::make_pair(left_exons[i].first,left_cds[0].first - 1)); }
                
                if (left_exons[i].first < left_cds[0].second && left_exons[i].second > left_cds[0].second )
                { utr3.push_back(std::make_pair(left_cds[0].second + 1,left_exons[i].second)); }
                if (left_exons[i].first > left_cds[0].second ) { utr3.push_back(left_exons[i]) ; }
            }
        }

        if (left_cds.size() == 2 ) {
            if (left_cds.size() == 2 ) {
                for(i=0 ; i < left_exons.size();i++) {
                    if (left_exons[i].second < left_cds[0].first) { utr5.push_back(left_exons[i]); }
                    if (left_exons[i].first < left_cds[0].first && left_exons[i].second > left_cds[0].first )
                    { utr5.push_back(std::make_pair(left_exons[i].first,left_cds[0].first - 1)); }
                    
                    if (left_exons[i].first < left_cds[1].second && left_exons[i].second > left_cds[1].second )
                    { utr3.push_back(std::make_pair(left_cds[1].second + 1,left_exons[i].second)); }
                    
                    if (left_exons[i].first > left_cds[1].second ) { utr5.push_back(left_exons[i]) ; }
                }
            }
        }
    }
    else {
        itgDn1k_st = std::max(long(0),exons[0].first - 1000);
        itgDn1k_end = std::max(long(0),exons[0].first -1);
        itgUp1k_st = exons[exons.size()-1].second + 1 ;
        itgUp1k_end = exons[exons.size()-1].second + 1000 ;
        itg1k.push_back(std::make_pair(itgUp1k_st,itgUp1k_end));
        itg1k.push_back(std::make_pair(itgDn1k_st,itgDn1k_end));

        if (left_cds.size() == 0 ) { utr5 = exons; };
        
        if (left_cds.size() == 1) {
            for( i=0;i < left_exons.size(); i++) {
                
                if (left_exons[i].second < left_cds[0].first) { utr3.push_back(left_exons[i]); }
                if (left_exons[i].first < left_cds[0].first && left_exons[i].second > left_cds[0].first ) { utr3.push_back(std::make_pair(left_exons[i].first,left_cds[0].first - 1)); }
                
                if (left_exons[i].first < left_cds[0].second && left_exons[i].second > left_cds[0].second )
                { utr5.push_back(std::make_pair(left_cds[0].second + 1,left_exons[i].second)); }
                if (left_exons[i].first > left_cds[0].second ) { utr5.push_back(left_exons[i]) ; }
            }
        }
        if (left_cds.size() == 2 ) {
            for(i=0 ; i < left_exons.size();i++) {
                if (left_exons[i].second < left_cds[0].first) { utr3.push_back(left_exons[i]); }
                if (left_exons[i].first < left_cds[0].first && left_exons[i].second > left_cds[0].first )
                { utr3.push_back(std::make_pair(left_exons[i].first,left_cds[0].first - 1)); }
                
                if (left_exons[i].first < left_cds[1].second && left_exons[i].second > left_cds[1].second )
                { utr5.push_back(std::make_pair(left_cds[1].second + 1,left_exons[i].second)); }

                if (left_exons[i].first > left_cds[1].second ) { utr5.push_back(left_exons[i]) ; }
            }
        }
        
    }

}



GeneFeatures::GeneFeatures(std::string GTFfilename,std::string id_attribute){

        read_features(GTFfilename,id_attribute);

}

GeneFeatures::~GeneFeatures(){
    chrom_itvTree_Dict_itr it;
    
    for (it=cds_exon_idx_plus.begin(); it != cds_exon_idx_plus.end(); it++) {
        std::map<int,IntervalTree *> tmp = it->second;
        std::map<int,IntervalTree *>::iterator tmp_itr ;
        //std::cout << tmp.size() << std::endl;
        for (tmp_itr = tmp.begin(); tmp_itr != tmp.end(); tmp_itr ++) {
            //std::cout << "delete tree " << std::endl;
            //std::cout << tmp_itr->first << std::endl;
            //std::cout << tmp_itr->second.size() << std::endl;
            delete tmp_itr->second;
        }
        
    }
    for (it=cds_exon_idx_minus.begin(); it != cds_exon_idx_minus.end(); it++) {
        std::map<int,IntervalTree *> tmp = it->second;
        std::map<int,IntervalTree *>::iterator tmp_itr ;
        
        for (tmp_itr = tmp.begin(); tmp_itr != tmp.end(); tmp_itr ++) {
            delete tmp_itr->second;
        }
        
    }
    
}

void GeneFeatures::build_tree(std::map<std::string, std::map<std::string,Gene> > temp_plus, std::map<std::string, std::map<std::string,Gene> > temp_minus)
{
    std::vector<Interval> itemlist;
    std::vector<Interval> sublist;
    std::map<std::string, std::map<std::string,Gene> >::iterator it;
    std::map<std::string,Gene> tmp;
    std::map<std::string,Gene>::iterator tmp_itr;
    gene_exon_Dict_It exon_str_itr;
    //Gene tmp_gene;
    
    std::string chr,gid;//,start_ss,end_ss;
    //std::stringstream st, end;
    //long st, end;
    int i, cur_bin_id,start_bin_id, js, je, k, buket_size;

    for (it = temp_plus.begin(); it != temp_plus.end(); it++) {
        chr = it->first;
        tmp = it->second;
        itemlist.clear();
        
        for (tmp_itr = tmp.begin(); tmp_itr != tmp.end(); tmp_itr++) {
            gid = tmp_itr->first;
            Gene tmp_gene = tmp_itr->second;
            tmp_gene.get_others();
            //std::cout << gid << " " << tmp_gene.exons.size() << std::endl;
            for (i=0; i< tmp_gene.exons.size(); i++) {
                std::string start_ss,end_ss;
                std::stringstream st, end;
                st << tmp_gene.exons[i].first;
                end << tmp_gene.exons[i].second;
                st >> start_ss;
                end >> end_ss;
                
                gene_exon_Dict_It gene_exon_it =plus_gene_exons.find(gid);
                
                if (gene_exon_it != plus_gene_exons.end()) {
                    std::vector<std::string> *tmp_ptr = &(plus_gene_exons[gid]);
                    (*tmp_ptr).push_back(start_ss+":"+end_ss+":+");
                }
                else {
                    std::vector<std::string> tmp_gene_str;
                    tmp_gene_str.push_back(start_ss+":"+end_ss+":+");
                    plus_gene_exons.insert(std::pair<std::string,std::vector<std::string> > (gid,tmp_gene_str));
                }
            }
            for (i=0; i < tmp_gene.cds.size(); i++) {
                itemlist.push_back(Interval(gid,tmp_gene.cds[i].first,tmp_gene.cds[i].second,"cds"));
                
            }
            for (i=0; i < tmp_gene.utr5.size(); i++) {
                itemlist.push_back(Interval(gid,tmp_gene.utr5[i].first,tmp_gene.utr5[i].second,"utr5"));
            }
            for (i=0; i < tmp_gene.utr3.size(); i++) {
                itemlist.push_back(Interval(gid,tmp_gene.utr3[i].first,tmp_gene.utr3[i].second,"utr3"));
            }
            for (i=0; i < tmp_gene.intron.size(); i++) {
                itemlist.push_back(Interval(gid,tmp_gene.intron[i].first,tmp_gene.intron[i].second,"intron"));
            }
            itemlist.push_back(Interval(gid,tmp_gene.itg1k[0].first,tmp_gene.itg1k[0].second,"Up1k"));
            itemlist.push_back(Interval(gid,tmp_gene.itg1k[1].first,tmp_gene.itg1k[1].second,"Dn1k")) ;
            
        }
        
        std::sort(itemlist.begin(),itemlist.end(),itv_comp) ; //key=operator.attrgetter('start'));
        //quick_sort(itemlist,0,itemlist.size()-1);
        start_bin_id = itemlist[0].start/BIN_SIZE;
        js = 0 ;
        je = 0 ;
        k = 0 ;
        
        
        for (i=0; i < itemlist.size(); i++) {
            cur_bin_id = itemlist[i].start/BIN_SIZE;
            //std::cout << cur_bin_id   << std::endl;
            
            if (cur_bin_id == start_bin_id) {
                je += 1;
            }
            else {
                buket_size = (int)sqrt(je - js) + 1;
                //std::cout << buket_size << std::endl;
                
                sublist = std::vector<Interval>(itemlist.begin()+js,itemlist.begin() + je);
                //std::cout << "sublist size " << sublist.size() << std::endl;
                /*if (start_bin_id == 0) {
                    for (int j=0;j < sublist.size(); j++) {
                        std::cout << "gene " << sublist[j].gene << " start " << sublist[j].start  << " end " << sublist[j].stop << " type " << sublist[j].type << std::endl;
                    }
                }*/
                //cds_exon_idx_plus[chr][start_bin_id] = new IntervalTree(sublist,16,buket_size,-1,-1,buket_size);
                cds_exon_idx_plus[chr][start_bin_id] = new IntervalTree(sublist);
               // std::cout << start_bin_id << " built one tree."  << std::endl;
                k ++;
                start_bin_id = cur_bin_id ;
                js = je ;
                je ++;
            }
        }

        
        if (js != je) {
            buket_size = (int) sqrt(je - js) + 1;
            sublist = std::vector<Interval>(itemlist.begin()+js,itemlist.begin() + je);
            //if (temp_gene.id == "DDX11L1") {
           // std::cout << buket_size << std::endl;
            
            
            //std::cout << "sublist size " << sublist.size() << std::endl;
            //}
            //cds_exon_idx_plus[chr][start_bin_id] = new IntervalTree(sublist, 16, buket_size,-1,-1,buket_size );
            cds_exon_idx_plus[chr][start_bin_id] = new IntervalTree(sublist);
        //print("tree depth = " + str(cds_exon_idx_plus[chr][start_bin_id].get_depth())+ "\n")
            k+=1;
        }
        
        
    }
    //std::cout << " build minus strand."  << std::endl;
    for (it = temp_minus.begin(); it != temp_minus.end(); it++) {
        //std::cout << "negative strand " << std::endl;
        chr = it->first;
        tmp = it->second;
        itemlist.clear();
        
        for (tmp_itr = tmp.begin(); tmp_itr != tmp.end(); tmp_itr++) {
            gid = tmp_itr->first;
            //std::cout << "GID  " << gid << std::endl;
            Gene tmp_gene = tmp_itr->second;
            tmp_gene.get_others();
            
            for (i=0; i< tmp_gene.exons.size(); i++) {
                std::string start_ss,end_ss;
                std::stringstream st, end;
                st << tmp_gene.exons[i].first;
                end << tmp_gene.exons[i].second;
                st >> start_ss;
                end >> end_ss;
                //minus_gene_exons[gid].push_back(start_ss+":"+end_ss+":-");
                
                gene_exon_Dict_It gene_exon_it =minus_gene_exons.find(gid);
                
                if (gene_exon_it != minus_gene_exons.end()) {
                    std::vector<std::string> *tmp_ptr = &(minus_gene_exons[gid]);
                    (*tmp_ptr).push_back(start_ss+":"+end_ss+":-");
                }
                else {
                    std::vector<std::string> tmp_gene_str;
                    tmp_gene_str.push_back(start_ss+":"+end_ss+":-");
                    minus_gene_exons.insert(std::pair<std::string,std::vector<std::string> > (gid,tmp_gene_str));
                }

                
            }
            
            for (i=0; i < tmp_gene.cds.size(); i++) {
                itemlist.push_back(Interval(gid,tmp_gene.cds[i].first,tmp_gene.cds[i].second,"cds"));
            }
            for (i=0; i < tmp_gene.utr5.size(); i++) {
                itemlist.push_back(Interval(gid,tmp_gene.utr5[i].first,tmp_gene.utr5[i].second,"utr5"));
            }
            for (i=0; i < tmp_gene.utr3.size(); i++) {
                itemlist.push_back(Interval(gid,tmp_gene.utr3[i].first,tmp_gene.utr3[i].second,"utr3"));
            }
            for (i=0; i < tmp_gene.intron.size(); i++) {
                itemlist.push_back(Interval(gid,tmp_gene.intron[i].first,tmp_gene.intron[i].second,"intron"));
            }
            itemlist.push_back(Interval(gid,tmp_gene.itg1k[0].first,tmp_gene.itg1k[0].second,"Up1k"));
            itemlist.push_back(Interval(gid,tmp_gene.itg1k[1].first,tmp_gene.itg1k[1].second,"Dn1k")) ;
            
        }
        
        //std::sort(key=operator.attrgetter('start'));
        //quick_sort(itemlist,0,itemlist.size()-1);
        std::sort(itemlist.begin(),itemlist.end(),itv_comp) ;
        start_bin_id = itemlist[0].start/BIN_SIZE;
        js = 0 ;
        je = 0 ;
        k = 0 ;
        
        for (i=0; i < itemlist.size(); i++) {
            cur_bin_id = itemlist[i].start/BIN_SIZE;
            if (cur_bin_id == start_bin_id) {
                je += 1;
            }
            else {
                buket_size = (int)sqrt(je - js) + 1;
               // std::cout << buket_size << std::endl;
                
                
                sublist = std::vector<Interval>(itemlist.begin()+js,itemlist.begin() + je);
                //cds_exon_idx_minus[chr][start_bin_id] = new IntervalTree(sublist,16,buket_size,-1,-1,buket_size);
                cds_exon_idx_minus[chr][start_bin_id] = new IntervalTree(sublist);
                k ++;
                start_bin_id = cur_bin_id ;
                js = je ;
                je ++;
            }
        }
        
        
        if (js != je) {
            buket_size = (int) sqrt(je - js) + 1;
            //std::cout << buket_size << std::endl;
            
            sublist = std::vector<Interval>(itemlist.begin()+js,itemlist.begin() + je);
            //cds_exon_idx_minus[chr][start_bin_id] = new IntervalTree(sublist, 16, buket_size,-1,-1,buket_size );
            cds_exon_idx_minus[chr][start_bin_id] = new IntervalTree(sublist);
            //print("tree depth = " + str(cds_exon_idx_plus[chr][start_bin_id].get_depth())+ "\n")
            k+=1;
        }
        
        
    }

}
//Reading & processing annotation files
void GeneFeatures::read_features(std::string gff_filename, std::string id_attribute)
{

    //dict of dicts since the builtin type doesn't support it for some reason
    std::map<std::string, std::map<std::string, Gene> > temp_plus ;
    std::map<std::string, std::map<std::string, Gene> > temp_minus ;
    std::map<std::string, std::map<std::string, Gene> >::iterator tmp_itr;
    std::map<std::string, Gene>::iterator id_itr;
    //std::map<std::string, Gene> gene_id_map;
    //std::vector<std::string> ::iterator id_itr;
    //Gene *g;
    //std::regex re("\\s*([^\\s\\=]+)[\\s=]+(.*)");
    //std::smatch sm;
    bool matched = false;
    //int k = 0;
    int i =  0;
    int counts = 0 ;
    int line_no = 0;
    long start = -1;
    long end = -1;
    std::size_t pos,cur_pos;
    //std::string left_str,sub_str;
    
    std::ifstream input; //(gff_filename);
    
    try{
        input.open (gff_filename, std::ifstream::in);
    
    //if (!std::getline(input,line)) {
    //    std::cout << "error read in file." << std::endl;
    //}
    while(! input.eof()){
    //while (std::getline(input,line)) {
        std::string line,chrom,source,feature,start_ss,end_ss,score,strand,frame,attributeStr;
        std::stringstream ss;
        std::string id = "";

        if (! std::getline(input,line)){
            break;
        }
        //std::cout << "read" << std::endl;
        line_no ++;
        
        if (line == "\n" || !line.compare(0,1,"#")) {
            continue;
        }
       // std::cout << line_no << std::endl;
        ss << line;
        std::getline(ss,chrom,'\t');
        
        std::getline(ss,source,'\t');
        std::getline(ss,feature,'\t');
        std::getline(ss,start_ss,'\t');
        std::getline(ss,end_ss,'\t');
        std::getline(ss,score,'\t');
        std::getline(ss,strand,'\t');
        std::getline(ss,frame,'\t');
        std::getline(ss,attributeStr,'\t');
        
        //std::cout << strand << std::endl;
        try{
        start = std::stol(start_ss);
        end = std::stol(end_ss);
        }
        catch (const std::invalid_argument& ia) {
            std::cerr << "Invalid argument: " << ia.what() << '\n';
            std::exit(1);
           // std::cout << start_ss << std::endl;
            //std::cout << end_ss << std::endl;
        }
        //( seqname, source, feature, start, end, score, strand, frame, attributeStr ) = line.split("\t")
        //id = self.__parse_GFF_attr_string(attributeStr,self.id_attribute);
        
        //pos = attributeStr.find(";");
        //std::string attr_items[] = attributeStr.split(";");
        cur_pos = 0;
        //std::cout <<attributeStr << std::endl;
        while (cur_pos < attributeStr.length()) {
            std::size_t next_pos = attributeStr.find(";",cur_pos);
            if (next_pos !=std::string::npos) {
                //std::cout << cur_pos << std::endl;
                //std::cout << next_pos << std::endl;
                
                std::string tok = attributeStr.substr(cur_pos,(next_pos - cur_pos));
                //std::cout << tok << std::endl;
                tok = tok.substr(tok.find_first_not_of(' '));
                pos = tok.find('=');
                if (pos == std::string::npos) {
                    pos = tok.find(' ');
                }
                std::string key = tok.substr(0,pos);
                key = key.substr(0,key.find(' '));
                
                std::size_t pos_stop = 1 + pos + tok.substr(pos+1).find_first_not_of(' ');
                std::string val = (tok[pos_stop] == '"') ? tok.substr(pos_stop+1, (tok.length() - (pos_stop+2))) : tok.substr(pos_stop, (tok.length() - (pos_stop+1)));
                
                if (key == id_attribute) {
                    id = val;
                    matched = true;
                    break;
                }
                cur_pos = next_pos + 1;
            }
            else {
                break;
            }
            
            
        }
        
        /*for (auto tok : tokens){
            //tok = attr_items[k].substr(attr_items[k].find_first_not_of(' '));
            tok = tok.substr(tok.find_first_not_of(''));
            pos = tok.find('=');
            if (pos == tok.npos) {
                pos = tok.find(' ');
            }
            std::string key = tok.substr(0,pos);
            key = key.substr(0,key.find(' '));
            
            std::size_t pos_stop = 1 + pos + tok.substr(pos+1).find_first_not_of(' ');
            std::string val = (tok[pos_stop] == '"') ? tok.substr(pos_stop+1, (tok.length() - (pos_stop+2))) : tok.substr(pos_stop, (tok.length() - (pos_stop+1)));
            
            if (key == id_attribute) {
                id = val;
                matched = true;
                break;
            }
        }*/
        /*matched = std::regex_match(attributeStr.substr(0,pos),sm,re);
        
        while (matched && pos + 1 != attributeStr.length() && id == "")
        {
            for (k=0; k < sm.size(); k++) {
                if (sm[k] == id_attribute) {
                    id = sm[k+1];
                    break;
                }
            }
            
            std::string left_str = attributeStr.substr(pos+1);
            cur_pos = pos;
            pos = left_str.find(";");
            matched = std::regex_match(attributeStr.substr(cur_pos,pos),sm,re);
            
        }*/
        if (!matched) {
        
            std::cout << "Failure parsing GFF attribute line." << std::endl;
            exit (EXIT_FAILURE);
        }
        //std::cout << "ID " << id << std::endl;
        
        if (id == "") {
            continue;
        }
        if ( std::find(features.begin(), features.end(), id) == features.end() ) {
            features.push_back(id);
        }
        
        if (feature == "CDS" ){
            if (strand == "+" ){
                tmp_itr = temp_plus.find(chrom);
                if (tmp_itr != temp_plus.end()) {
                    id_itr = temp_plus[chrom].find(id);
                    if (id_itr != temp_plus[chrom].end()) {
                        Gene *g = &(id_itr->second);
                        //(id_itr->second).add_cds(start,end);
                        (*g).add_cds(start,end);
                    }
                    else{
                        Gene g (id,strand);
                        g.add_cds(start,end);
                        temp_plus[chrom].insert(std::pair<std::string,Gene>(id,g));
                        //std::map<std::string, Gene>  gene_id_map ;
                        //gene_id_map.insert(std::pair<std::string,Gene> (id,g));
                        //temp_plus[chrom] = gene_id_map;
                        //temp_plus[chrom].insert(gene_id_map);
                        //temp_plus.insert(std::pair<std::string,std::map<std::string, Gene> > (chrom,gene_id_map));
                    }
                }
                else{
                    Gene g(id,strand);
                    g.add_cds(start,end);
                    std::map<std::string, Gene>  gene_id_map ;
                    gene_id_map.insert(std::pair<std::string,Gene> (id,g));
                    //gene_id_map[id] = g;
                    //temp_plus[chrom] = gene_id_map;
                    //temp_plus[chrom].insert(gene_id_map);
                    temp_plus.insert(std::pair<std::string,std::map<std::string, Gene> > (chrom,gene_id_map));
                    
                }

            }
                    
            if (strand == "-" ) {
                
                tmp_itr = temp_minus.find(chrom);
                if (tmp_itr != temp_minus.end()) {
                    id_itr = temp_minus[chrom].find(id);
                    if (id_itr != temp_minus[chrom].end()) {
                        Gene *g = &(id_itr->second);
                        (*g).add_cds(start,end);
                    }
                    else{
                        Gene g(id,strand);
                        g.add_cds(start,end);
                        std::map<std::string, Gene>  gene_id_map ;
                        temp_minus[chrom].insert(std::pair<std::string,Gene>(id,g));
                        //gene_id_map.insert(std::pair<std::string,Gene> (id,g));
                        //gene_id_map[id] = g;
                        //temp_plus[chrom] = gene_id_map;
                        //temp_minus[chrom].insert(gene_id_map);
                        //temp_minus.insert(std::pair<std::string,std::map<std::string, Gene> > (chrom,gene_id_map));
                        //gene_id_map[id] = g;
                        //temp_minus[chrom] = gene_id_map;
                        //temp_minus[chrom].insert(std::map<std::string,Gene>(id,g));
                    }
                }
                else{
                    Gene g(id,strand);
                    g.add_cds(start,end);
                    std::map<std::string, Gene> gene_id_map ;
                    gene_id_map.insert(std::pair<std::string,Gene> (id,g));
                    //gene_id_map[id] = g;
                    //temp_minus[chrom] = gene_id_map;
                    temp_minus.insert(std::pair<std::string,std::map<std::string, Gene> > (chrom,gene_id_map));
                    //temp_minus[chrom].insert(std::map<std::string,Gene>(id,g));
                    
                    
                }

            }
        }
        
        if (feature == "exon" ){
            counts += 1 ;
            if (strand == "+" ){
                tmp_itr = temp_plus.find(chrom);
                if (tmp_itr != temp_plus.end()) {
                    id_itr = temp_plus[chrom].find(id);
                    if (id_itr != temp_plus[chrom].end()) {
                        Gene *g = &(id_itr->second);
                        (*g).add_exons(start,end);
                    }
                    else{
                        Gene g (id,strand);
                        g.add_exons(start,end);
                        //std::map<std::string, Gene> gene_id_map ;
                        //gene_id_map[id] = g;
                        //temp_plus[chrom] = gene_id_map;
                        //temp_plus[chrom].insert(std::map<std::string,Gene>(id,g));
                        //gene_id_map.insert(std::pair<std::string,Gene> (id,g));
                        //temp_plus.insert(std::pair<std::string,std::map<std::string, Gene> > (chrom,gene_id_map));
                        temp_plus[chrom].insert(std::pair<std::string,Gene>(id,g));
                    }
                }
                else{
                    Gene g (id,strand);
                    g.add_exons(start,end);
                    std::map<std::string, Gene> gene_id_map ;
                    //gene_id_map[id] = g;
                    //temp_plus[chrom] = gene_id_map;
                    gene_id_map.insert(std::pair<std::string,Gene> (id,g));
                    temp_plus.insert(std::pair<std::string,std::map<std::string, Gene> > (chrom,gene_id_map));
                    //temp_plus[chrom].insert(std::map<std::string,Gene>(id,g));
                    
                }
                
            }
            
            if (strand == "-" ) {
               // std::cout << "negative strand" << std::endl;
                //std::cout << id << std::endl;
                tmp_itr = temp_minus.find(chrom);
                if (tmp_itr != temp_minus.end()) {
                    
                    id_itr = temp_minus[chrom].find(id);
                    if (id_itr != temp_minus[chrom].end()) {
                        Gene *g = &(id_itr->second);
                        (*g).add_exons(start,end);
                    }
                    else{
                        Gene g(id,strand);
                        g.add_exons(start,end);
                        //std::map<std::string, Gene> gene_id_map ;
                        //gene_id_map.insert(std::pair<std::string,Gene> (id,g));
                        //gene_id_map[id] = g;
                        //temp_minus[chrom] = gene_id_map;
                        temp_minus[chrom].insert(std::pair<std::string,Gene>(id,g));
                        //temp_minus.insert(std::pair<std::string,std::map<std::string, Gene> > (chrom,gene_id_map));
                        //temp_minus[chrom]
                    }
                }
                else{
                    Gene g (id,strand);
                    g.add_exons(start,end);
                    std::map<std::string, Gene> gene_id_map ;
                    gene_id_map.insert(std::pair<std::string,Gene> (id,g));
                    //gene_id_map[id] = g;
                    //temp_minus[chrom] = gene_id_map;
                    //temp_minus[chrom].insert(std::map<std::string,Gene>(id,g));
                    //std::cout << "insert to negative strand" << std::endl;
                    temp_minus.insert(std::pair<std::string,std::map<std::string, Gene> > (chrom,gene_id_map));
                    
                }
                
            }
            
                    
        }
    
        i += 1 ;
        if (i % 100000 == 0 )
        {
            //sys.stderr.write("%d GTF lines processed.\n" % i);
            //std::cout << i << " GTF lines processed." << std::endl;
        }
        
    }
    
    
    input.close();
    //    except:
      //      sys.stderr.write("Error occured in %s.\n" % gff.get_line_number_string())
        //raise
        
    }
    catch(std::ifstream::failure e){
        std::cout << "error in read file " << gff_filename << std::endl;
    }

    if (counts == 0 ){
        std::cout << "Warning: No features of type 'exon' or 'CDS' found in gene GTF file." << std::endl;
    }


        //build interval trees
        //trees for CDS
        //sys.stderr.write("start to build interval tree\n")
    //std::cout << "temp_plus size " << temp_plus.size() << std::endl;
    //std::cout << "temp_plus chrom1 size" << temp_plus["chr1"].size() << std::endl;
    build_tree(temp_plus,temp_minus);
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
std::vector<std::pair<long,long> > GeneFeatures::find_exons(std::string chrom,long st,long end,std::string strand,std::string gene)
{
    std::vector<std::pair<long,long> > exons;
    chrom_itvTree_Dict::iterator chrom_it;
    std::map<int,IntervalTree *>::iterator bin_iter;
    std::vector<Interval> fs ;
    std::vector<Interval> temp ;
    //std::vector<Interval> fs1 ;
    int i, bin_id ;
    
    //try:
    bin_id = st/BIN_SIZE ;
    
    if (strand == "+" || strand == "."){
        chrom_it = cds_exon_idx_plus.find(chrom);
        if (chrom_it != cds_exon_idx_plus.end()) {
            bin_iter = cds_exon_idx_plus[chrom].find(bin_id);
            if (bin_iter != cds_exon_idx_plus[chrom].end()) {
                fs = (*cds_exon_idx_plus[chrom][bin_id]).find(st,end);
            }
        }
    }
            
    if (strand == "-" || strand == "."){
        chrom_it = cds_exon_idx_minus.find(chrom);
        if (chrom_it != cds_exon_idx_minus.end()) {
            bin_iter = cds_exon_idx_minus[chrom].find(bin_id);
            if (bin_iter != cds_exon_idx_minus[chrom].end()) {
                temp = (*cds_exon_idx_minus[chrom][bin_id]).find(st,end);
                fs.insert(fs.end(),temp.begin(),temp.end());
            }
        }
    }
    
    //    except:
      //      raise
        
    for(i =0 ; i < fs.size(); i++){
        if (fs[i].gene == gene && fs[i].type != "intron" && fs[i].type != "Up1k" && fs[i].type != "Dn1k"){
            exons.push_back(std::make_pair(fs[i].start,fs[i].stop));
        }
    }

    std::sort(exons.begin(),exons.end());
    return exons;
}

std::vector<std::string> GeneFeatures::getFeatures() {
    return features ;
}

std::vector<std::string> GeneFeatures::get_exon_pair(std::string g){
    std::vector<std::string> res;
    gene_exon_Dict_It it = plus_gene_exons.find(g);
    
    if (it != plus_gene_exons.end() ) {
        return plus_gene_exons[g];
    }
    else{
        it =minus_gene_exons.find(g);
        if (it != minus_gene_exons.end()) {
            return minus_gene_exons[g];
        }
        else {
            return res;
        }
    }
}
    
std::vector<std::string> GeneFeatures::Gene_annotation(std::string chrom, std::vector<std::pair<long,long> > itv_list,std::string strand)
{
    std::vector<std::string> genes;
    std::vector<Interval> fs ;
   // std::vector<Interval> tmp ;
    chrom_itvTree_Dict::iterator chrom_it;
    std::map<int,IntervalTree *>::iterator bin_iter;
    
    std::string type="" ;
    int i,bin_id_s, bin_id_e;
    bool isCDS = false;
    bool isINTRON = false;
    bool isUTR3 = false;
    bool isUTR5 = false;
    bool isITGup = false;
    bool isITGdn = false;
    
    for (i=0; i < itv_list.size(); i ++) {
        bin_id_s = itv_list[i].first/BIN_SIZE;
        bin_id_e = itv_list[i].second/BIN_SIZE;
        
        if (strand == "+" || strand == ".") {
            chrom_it = cds_exon_idx_plus.find(chrom);
            if (chrom_it !=  cds_exon_idx_plus.end()) {
                
                bin_iter = cds_exon_idx_plus[chrom].find(bin_id_s);
                if (bin_iter != cds_exon_idx_plus[chrom].end()  )
                {
                    //std::cout << "start to search tree" << std::endl;
                    fs = (*cds_exon_idx_plus[chrom][bin_id_s]).find(itv_list[i].first,itv_list[i].second);
                    if (bin_id_s != bin_id_e) {
                        bin_iter = cds_exon_idx_plus[chrom].find(bin_id_e);
                        if (bin_iter != cds_exon_idx_plus[chrom].end()  ){
                        std::vector<Interval>  tmp = (*cds_exon_idx_plus[chrom][bin_id_e]).find(itv_list[i].first,itv_list[i].second);
                            //for (int j =0; j < tmp.size(); j++) {
                          //      fs.push_back(tmp[j]);
                        //    }
                        fs.insert(fs.end(),tmp.begin(),tmp.end());
                        }
                    }
                    //std::cout << fs.size() << std::endl;
                }
            }
        }

        if (strand == "-" or strand == ".") {
            chrom_it = cds_exon_idx_minus.find(chrom);
            if (chrom_it !=  cds_exon_idx_minus.end()) {
                
                bin_iter = cds_exon_idx_minus[chrom].find(bin_id_s);
                if (bin_iter != cds_exon_idx_minus[chrom].end()  )
                {
                    //std::cout << "start to search tree" << std::endl;
                    std::vector<Interval>  tmp = (*cds_exon_idx_minus[chrom][bin_id_s]).find(itv_list[i].first,itv_list[i].second);
                    //for (int j =0; j < tmp.size(); j++) {
                      //  fs.push_back(tmp[j]);
                    //}
                    fs.insert(fs.end(),tmp.begin(),tmp.end());

                    if (bin_id_s != bin_id_e) {
                        bin_iter = cds_exon_idx_minus[chrom].find(bin_id_e);
                        if (bin_iter != cds_exon_idx_minus[chrom].end()  ){
                        std::vector<Interval>  tmp = (*cds_exon_idx_minus[chrom][bin_id_e]).find(itv_list[i].first,itv_list[i].second);
                          //  for (int j =0; j < tmp.size(); j++) {
                            //    fs.push_back(tmp[j]);
                            //}
                        fs.insert(fs.end(),tmp.begin(),tmp.end());
                        }
                    }
                    //std::cout << fs.size() << std::endl;
                    
                }
           }
        }
    }
    //std::cout << fs.size() << std::endl;
    for(i =0 ; i < fs.size(); i++){
        //std::cout << fs[i].type << "\t" << fs[i].gene << std::endl;
        if (fs[i].type == "cds") {
            isCDS = true ;
            //type = "cds";
            break;
        }
        if (fs[i].type == "utr5") {
            isUTR5 = true;
            
        }
        if (fs[i].type == "utr3") {
            isUTR3 = true;
            
        }
        if (fs[i].type == "intron") {
            isINTRON = true;
            
        }
        if (fs[i].type == "Up1k") {
            isITGup = true;
            
        }
        if (fs[i].type == "Dn1k") {
            isITGdn = true;
        }
        
    }
    
    if (isCDS) {
        type = "cds";
    }
    else{
     if (isUTR5){
        type = "utr5";
     }
    else { if (isUTR3){
        type = "utr3";
      }
    else { if (isINTRON) {
        type = "intron";
      }
    else { if (isITGup) {
        type = "Up1k";
      }
    else { if (isITGdn) {
        type = "Dn1k";
      }
    else {
        return genes;
     }
    }}}}}
     
    
    genes.push_back(type);
    for (i=0; i < fs.size(); i++) {
        if (fs[i].type == type) {
            genes.push_back(fs[i].gene);
        }
    }
    return genes;

}

/*
int main() {
    std::string filename = "./test.gtf";
    std::string id_attr = "gene_id";
   
    std::vector<chr_ITV> itv_list ;
    chr_ITV exp ;
    exp.chrom = "chr1";
    exp.start = 11870;
    exp.end = 73220;
    itv_list.push_back(exp);
    
    std::cout << "start to build tree " << std::endl;
    GeneFeatures gIdx (filename,id_attr);
    std::cout << "after  build tree " << std::endl;

   // for (int i=0; i < itv_list.size(); i++) {
    //    std::cout << itv_list[i].start << std::endl;
    //}
    std::vector<std::string> res = gIdx.Gene_annotation(itv_list,".");
    for (int i=0;i<res.size(); i++) {
        std::cout << res[i] << std::endl;
    }
    
}*/
