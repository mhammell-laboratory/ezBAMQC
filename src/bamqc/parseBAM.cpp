//
//  parseBAM.cpp
//  BAMQC_c++
//
//  Created by Ying Jin on 11/11/15.
//  Copyright (c) 2015 Ying Jin. All rights reserved.
//

#include "parseBAM.h"
#include "Coverage_prof.h"
#include "GeneFeatures.h"
#include "InnerDist_prof.h"
#include "Mappability.h"
#include "Results.h"
#include "rRNA.h"
#include "Constants.h"

//#include <malloc.h>
#include <string>
#include <stdlib.h>
#include <pthread.h>
//#include <stdio.h>
#include <iostream>
#include <fstream>
#include <unistd.h>

#include "htslib/sam.h"



unsigned int tick_time =3000;

#define IS_PAIRED(bam) ((bam)->core.flag&BAM_FPAIRED)
#define IS_QCFAIL(bam) ((bam)->core.flag & BAM_FQCFAIL)

#define IS_PAIRED_AND_MAPPED(bam) (((bam)->core.flag&BAM_FPAIRED) && !((bam)->core.flag&BAM_FUNMAP) && !((bam)->core.flag&BAM_FMUNMAP))

#define IS_PROPERLYPAIRED(bam) (((bam)->core.flag&(BAM_FPAIRED|BAM_FPROPER_PAIR)) == (BAM_FPAIRED|BAM_FPROPER_PAIR) && !((bam)->core.flag&BAM_FUNMAP))

#define IS_UNMAPPED(bam) ((bam)->core.flag&BAM_FUNMAP)
#define IS_REVERSE(bam) ((bam)->core.flag&BAM_FREVERSE)
#define IS_MATE_REVERSE(bam) ((bam)->core.flag&BAM_FMREVERSE)
#define IS_READ1(bam) ((bam)->core.flag&BAM_FREAD1)
#define IS_READ2(bam) ((bam)->core.flag&BAM_FREAD2)
#define IS_DUP(bam) ((bam)->core.flag&BAM_FDUP)
//#defind READ_NAME(bam) (bam_get_qname(bam))

/*DEF BAM_CMATCH     = 0
DEF BAM_CINS       = 1
DEF BAM_CDEL       = 2
DEF BAM_CREF_SKIP  = 3
DEF BAM_CSOFT_CLIP = 4
DEF BAM_CHARD_CLIP = 5
DEF BAM_CPAD       = 6
DEF BAM_CEQUAL     = 7
DEF BAM_CDIFF      = 8*/

typedef struct
{
    bam1_t * first;
    bam1_t * second;
    int type;
} read_pair_t ;

typedef struct
{
    unsigned short thread_id;
    pthread_t thread_object;
    //pthread_spinlock_t cur_reads_lock;
    pthread_spinlock_t cur_reads_lock;
    
    //std::vector<std::pair<bam1_t *,bam1_t *> > *cur_reads;
    std::vector<read_pair_t> * cur_reads;
    
    Results * res;
    InnerDist_prof * inDist_prof;
    Clipping_prof * clip_prof;
    Coverage_prof * cov_prof;
    
} thread_context_t;

typedef struct
{
    unsigned short thread_number;
    int all_finished;
    std::string format;
    std::string stranded;
    //bam_hdr_t * header;
	std::vector<std::string> refnames;
    std::vector<thread_context_t *> thread_contexts;
    
    GeneFeatures * geneIdx;
    rRNA * rRNAIdx;
    
} global_context_t;

struct arg_struct {
    global_context_t * arg1;
	thread_context_t * arg2;
};
std::vector<std::pair<int, int> > fetch_intron(int st, uint32_t * cigar, uint32_t n_cigar,std::string format)
{
    //''' fetch intron regions defined by cigar. st must be zero based return list of tuple of (st, end)'''
    //match = re.compile(r'(\d+)(\D)')
    int chrom_st = st;
    if (format == "BAM") { chrom_st += 1 ;}
    
    std::vector<std::pair<int,int> > intron_bound ;
    
    for (unsigned int i=0; i < n_cigar; i++){ //code and size
        if (bam_cigar_op(cigar[i])==BAM_CMATCH) { chrom_st += bam_cigar_oplen(cigar[i]);} //match
        if (bam_cigar_op(cigar[i])==BAM_CINS) {continue;} //insertion
        if (bam_cigar_op(cigar[i])==BAM_CDEL) { chrom_st += bam_cigar_oplen(cigar[i]);} //deletion
        if (bam_cigar_op(cigar[i])==BAM_CREF_SKIP) { intron_bound.push_back(std::pair<int,int> (chrom_st,chrom_st+ bam_cigar_oplen(cigar[i])-1)); } //gap
        if (bam_cigar_op(cigar[i])==BAM_CSOFT_CLIP) { chrom_st += bam_cigar_oplen(cigar[i]);}// soft clipping
    }
    
    return intron_bound ;
}

std::vector<std::pair<int, int> > fetch_exon(int st,uint32_t * cigar,uint32_t n_cigar,std::string format)
{
    //''' fetch exon regions defined by cigar. st must be zero based return list of tuple of (st, end)'''
    //match = re.compile(r'(\d+)(\D)')
    int chrom_st = st;
    if (format == "BAM") { chrom_st += 1;}
    
    std::vector<std::pair<int,int> >exon_bound;

    for (unsigned int i=0; i < n_cigar; i++){ //code and size
        if (bam_cigar_op(cigar[i])==BAM_CMATCH) { exon_bound.push_back(std::pair<int, int> (chrom_st,chrom_st + bam_cigar_oplen(cigar[i])-1));} //match
        if (bam_cigar_op(cigar[i])==BAM_CINS) {continue;} //insertion
        if (bam_cigar_op(cigar[i])==BAM_CDEL) { chrom_st += bam_cigar_oplen(cigar[i]);} //deletion
        if (bam_cigar_op(cigar[i])==BAM_CREF_SKIP) { chrom_st += bam_cigar_oplen(cigar[i]); } //gap
        if (bam_cigar_op(cigar[i])==BAM_CSOFT_CLIP) { chrom_st += bam_cigar_oplen(cigar[i]);}// soft clipping
    }
    return exon_bound;
}


std::pair<int, int> ovp_gene( GeneFeatures * geneIdx,std::string chrom1,std::vector<std::pair<int,int> >exon_blocks1,std::string chrom2,std::vector<std::pair<int, int> > exon_blocks2,std::string strand,std::vector<int> * mapped_exons)
{
    std::map<int,int> res1 = geneIdx->Gene_annotation(chrom1,exon_blocks1,strand,mapped_exons);
    std::map<int,int> res2 = geneIdx->Gene_annotation(chrom2,exon_blocks2,strand,mapped_exons);
    if(res2.size() == 0) 
    {
        if(res1.size() ==0 ) { return std::pair<int,int>(-1,-1);}
        if(res1.size() == 1) { return std::pair<int,int>(res1.begin()->second,res1.begin()->first);}
        std::map<int,std::vector<int> > type_gene_list;
        int min_type = 6;
        int second_min_type = 6;
        for(auto& kv : res1){
             if(type_gene_list.find(kv.second) != type_gene_list.end())
             {
                  type_gene_list[kv.second].push_back(kv.first);
             }
             else {
                  std::vector<int> tmp;
                  tmp.push_back(kv.first);
                  type_gene_list.insert(std::pair<int,std::vector<int> >(kv.second,tmp));
             }
             if(kv.second < min_type ) {min_type = kv.second; second_min_type = min_type;}
        }
        if(type_gene_list[min_type].size() == 1 && min_type <= 3 && second_min_type > 3)
        {  return std::pair<int,int> (min_type,type_gene_list[min_type][0]); }
        return std::pair<int,int>(min_type,-1);
    }
    if(res1.size() == 0) 
    {
        if(res2.size() ==0 ) { return std::pair<int,int>(-1,-1);}
        if(res2.size() == 1) { return std::pair<int,int>(res2.begin()->second,res2.begin()->first);}
        std::map<int,std::vector<int> > type_gene_list;
        int min_type = 6;
        int second_min_type = 6;
        for(auto& kv : res2){
             if(type_gene_list.find(kv.second) != type_gene_list.end())
             {
                  type_gene_list[kv.second].push_back(kv.first);
             }
             else {
                  std::vector<int> tmp;
                  tmp.push_back(kv.first);
                  type_gene_list.insert(std::pair<int,std::vector<int> >(kv.second,tmp));
             }
             if(kv.second < min_type ) {min_type = kv.second; second_min_type = min_type;}
        }
        if(type_gene_list[min_type].size() == 1 && min_type <= 3 && second_min_type > 3)
        {  return std::pair<int,int> (min_type,type_gene_list[min_type][0]); }
        return std::pair<int,int>(min_type,-1);
    }
    //PE
    if (res1.size() ==1 && res2.size() ==1 ){
      int g1 = res1.begin()->first;
      int type1 = res1.begin()->second;
      int g2 = res2.begin()->first;
      int type2 = res2.begin()->second;
      if(g1 == g2) {
            if(type1 < type2) { return std::pair<int,int> (type1,g1); }
            else { return std::pair<int,int> (type2,g1); }
      }else {
            if(type1 < type2) { return std::pair<int,int> (type1,-1); }
            else { return std::pair<int,int> (type2,-1); }
      }
    }
    if (res1.size() > 1 || res2.size() > 1 )
    {
        std::vector<int> ovp_genes ;
        
        for(auto& kv : res1 )
        {
            if(res2.find(kv.first)!=res2.end()) {
                ovp_genes.push_back(kv.first);
            }
        }
        if(ovp_genes.size() == 1) {
            if(res1[ovp_genes[0]] < res2[ovp_genes[0]]) { return std::pair<int,int> (res1[ovp_genes[0]],ovp_genes[0]); }
            else { return std::pair<int,int> (res2[ovp_genes[0]],ovp_genes[0]); }
        }
        if (ovp_genes.size() > 1) {
            std::map<int,std::vector<int> > type_gene_list;
            int min_type = 6;
            int second_min_type = 6;
            for(auto g : ovp_genes){
                int type1 = res1[g];
                int type2 = res2[g];
                int type = type1;
                if (type2 < type1) {
                    type = type2;
                }
                if(type_gene_list.find(type) != type_gene_list.end())
                {
                    type_gene_list[type].push_back(g);
                }
                else {
                    std::vector<int> tmp;
                    tmp.push_back(g);
                    type_gene_list.insert(std::pair<int,std::vector<int> >(type,tmp));
                }
                if(type < min_type ) {min_type = type; second_min_type = min_type;}
            }
            if(type_gene_list[min_type].size() == 1 && min_type <= 3 && second_min_type > 3)
            {  return std::pair<int,int> (min_type,type_gene_list[min_type][0]); }
            else { return std::pair<int,int>(min_type,-1); }

        }
        
        if (ovp_genes.size() == 0) {
            std::map<int,std::vector<int> > type_gene_list;
            int min_type = 6;
            int second_min_type = 6;
            for (auto& kv : res1) {
                if(type_gene_list.find(kv.second) != type_gene_list.end())
                {
                    type_gene_list[kv.second].push_back(kv.first);
                }
                else {
                    std::vector<int> tmp;
                    tmp.push_back(kv.first);
                    type_gene_list.insert(std::pair<int,std::vector<int> >(kv.second,tmp));
                }
                if(kv.second < min_type ) {min_type = kv.second; second_min_type = min_type;}
            }
            for(auto& kv : res2){
                if(type_gene_list.find(kv.second) != type_gene_list.end())
                {
                    type_gene_list[kv.second].push_back(kv.first);
                }
                else {
                    std::vector<int> tmp;
                    tmp.push_back(kv.first);
                    type_gene_list.insert(std::pair<int,std::vector<int> >(kv.second,tmp));
                }
                if(kv.second < min_type ) {min_type = kv.second; second_min_type = min_type;}
            }
            if(type_gene_list[min_type].size() == 1 && min_type <= 3 && second_min_type > 3)
            {  return std::pair<int,int> (min_type,type_gene_list[min_type][0]); }
            else { return std::pair<int,int>(min_type,-1); }
            
            
        }
        
    }
    
    return std::pair<int, int> (-1,-1);
    
}

//void process_aligned_fragment(global_context_t * global_context,thread_context_t * thread_context,std::pair<bam1_t *,bam1_t *> read_pair)
void process_aligned_fragment(global_context_t * global_context,thread_context_t * thread_context,read_pair_t  read_pair)
{
    
    bam1_t * cur_read1 = read_pair.first;
    bam1_t * cur_read2 = read_pair.second;
    std::string strand1 = ".";
    std::string strand2 = ".";
    int chrom1_id = -1;
    int chrom2_id = -1;
    std::vector<std::pair<int, int> > exons1 ;
    std::vector<std::pair<int, int> > exons2 ;
    std::vector<std::pair<int, int> > intron_blocks1 ;
    std::vector<std::pair<int, int> > intron_blocks2 ;
    std::string qname = "";
	std::string chrom1 = "";
	std::string chrom2 = "";

    if(cur_read1 != NULL ) {
        if(!IS_UNMAPPED(cur_read1))
        {
            uint32_t *cigar = bam_get_cigar(cur_read1);
            char * name = bam_get_qname(cur_read1);
            
            if (read_pair.type == 0) { //only for uniq-read
                thread_context-> clip_prof->set(cur_read1->core.l_qseq,cur_read1->core.n_cigar,cigar,cur_read1->core.qual);
            }
            chrom1_id  =  cur_read1->core.tid;

            qname = std::string(name);
			if((size_t)cur_read1->core.tid < global_context->refnames.size()){
			chrom1 = global_context->refnames[cur_read1->core.tid];
			}
			else {
			  std::cout << (size_t)cur_read1->core.tid  << "\t" << "missing reference sequences." << std::endl;
			  std::exit(1);
			}

            strand1 = IS_REVERSE(cur_read1) ? "-" : "+";
            
            if (global_context->stranded =="reverse") {
                strand1 = (strand1 == "+") ? "-" : "+";
            }
            exons1 = fetch_exon(cur_read1->core.pos, cigar,cur_read1->core.n_cigar, global_context->format);
			intron_blocks1 = fetch_intron(cur_read1->core.pos,cigar,cur_read1->core.n_cigar,global_context->format);
        }
    }
    if(cur_read2 != NULL){
        if(! IS_UNMAPPED(cur_read2)){
            uint32_t *cigar = bam_get_cigar(cur_read2);
            if (read_pair.type ==0) {
                thread_context->clip_prof->set(cur_read2->core.l_qseq,cur_read2->core.n_cigar,cigar,cur_read2->core.qual);
            }
            strand2 =  IS_REVERSE(cur_read2) ? "-" : "+";
            
            if (global_context-> stranded =="reverse") {
                strand2 = ( strand2 == "+") ? "-" : "+";
            }
            chrom2_id = cur_read2->core.tid;
			if((size_t)cur_read2->core.tid < global_context->refnames.size()){
			chrom2 = global_context->refnames[cur_read2->core.tid];
			}
			else {
			  std::cout << (size_t)cur_read2->core.tid << "\t" <<"missing reference sequences." << std::endl;
			  std::exit(1);
			}
            exons2 = fetch_exon(cur_read2->core.pos, cigar,cur_read2->core.n_cigar, global_context->format);
            intron_blocks2 = fetch_intron(cur_read2->core.pos,cigar,cur_read2->core.n_cigar,global_context->format);
        }
    }

    if (cur_read1 == NULL || IS_UNMAPPED(cur_read1)) {
        strand1 =  (strand2 == "-") ? "+" : "-";
    }
    std::string strand = strand1;
    if (global_context->stranded == "no"){
        strand = ".";
    }

    //type =1 means multi-read. rRNA read
    if (read_pair.type == 1) {
        
        //std::cout << chrom1 << "\t" << chrom2 << std::endl;

        if (global_context->rRNAIdx != NULL){
            if (global_context->rRNAIdx->is_rRNA(chrom1,exons1,strand) || global_context->rRNAIdx->is_rRNA(chrom2,exons2,strand)){
                thread_context->res->rRNA_read += 1;
                //std::cout << qname << std::endl;
                if(cur_read1 != NULL){
                    bam_destroy1(cur_read1);
                }
                if(cur_read2 != NULL){
                    bam_destroy1(cur_read2);
                }
            return;
            }
        }
    }
    else { //uniq read

    
    if (intron_blocks1.size() + intron_blocks2.size() == 0){
        thread_context->res->noSplice += 1;
    }
    else {
        thread_context->res->splice += 1;
    }
    //paired read
    if (chrom1_id != -1 && chrom2_id != -1){
        thread_context->res->paired_reads += 1;
        
        if (strand1 == "-" && strand2 == "-" ){
            thread_context->res->mapped_minus_minus  += 1;
        }
        if (strand1 == "+" && strand2 == "+" ){
            thread_context->res->mapped_plus_plus += 1;
        }
        if (strand1 == "+" && strand2 == "-" ){
            thread_context->res->mapped_plus_minus += 1;
        }
        if (strand1 == "-" && strand2 == "+"){
            thread_context->res->mapped_minus_plus += 1;
        }
        if  (chrom1_id != chrom2_id ){
            thread_context->res->paired_diff_chrom += 1;
	        if(cur_read1 != NULL){
	            bam_destroy1(cur_read1);
	        }
	        if(cur_read2 != NULL){
                bam_destroy1(cur_read2);
	        }
            return;
        }
    }

    //mapped read len
    //rRNA read
    if (global_context->rRNAIdx != NULL){
        if (global_context->rRNAIdx->is_rRNA(chrom1,exons1,strand) || global_context->rRNAIdx->is_rRNA(chrom2,exons2,strand)){
            thread_context->res->rRNA_read += 1;
            //std::cout << qname << std::endl;
	    if(cur_read1 != NULL){
	       bam_destroy1(cur_read1);
	    }
	    if(cur_read2 != NULL){
	       bam_destroy1(cur_read2);
	    }
            return;
        }
    }
    //cur_time = time.time()
    std::vector<int> * mapped_exons = new std::vector<int>();
    
    std::pair<int,int> ovp = ovp_gene(global_context->geneIdx,chrom1,exons1,chrom2,exons2,strand,mapped_exons);
    switch(ovp.first)
    {
        case CDS:
            thread_context->res->cds_exon_read +=1;
            break;
        case UTR5:
            thread_context->res->utr_5_read +=1;
            break;
        case UTR3:
            thread_context->res->utr_3_read +=1;
            break;
        case INTRON:
            thread_context->res->intron_read +=1;
            ovp.second = -1;
            break;
        case ITGUP1K:
            thread_context->res->intergenic_up1kb_read +=1;
            ovp.second = -1;
            break;
        case ITGDN1K:
            thread_context->res->intergenic_down1kb_read +=1;
            ovp.second = -1;
            break;
        default:
            thread_context->res->intergenic_read +=1;
            break;
    }
    thread_context->cov_prof->count(ovp.second,exons1,exons2,*mapped_exons);
    if(chrom1 != "" && chrom2 != ""){
        if( IS_PROPERLYPAIRED(cur_read1)){
            if (strand1 == "+" || strand1 == "."){
                if (global_context->stranded == "reverse") {
                    thread_context->inDist_prof->count(global_context->geneIdx,ovp.first,cur_read2,chrom2,intron_blocks2,strand1);
                }
                else{
                    thread_context->inDist_prof->count(global_context->geneIdx,ovp.first,cur_read1,chrom1,intron_blocks1,strand1);
                }
            }
            else {
                if (global_context->stranded == "reverse" ){
                    thread_context->inDist_prof->count(global_context->geneIdx,ovp.first,cur_read1,chrom1,intron_blocks1,strand1);
                }
                else {
                    thread_context->inDist_prof->count(global_context->geneIdx,ovp.first,cur_read2,chrom2,intron_blocks2,strand1);
                }
            }
        }
    }
    delete mapped_exons;
	if(cur_read1 != NULL){
	    bam_destroy1(cur_read1);
	}
	if(cur_read2 != NULL){
	    bam_destroy1(cur_read2);
	}
  }//uniq read
}

void* worker(void * vargs)
{
    struct arg_struct * args = (struct arg_struct *) vargs;
    thread_context_t * thread_context = args->arg2;
    global_context_t * global_context = args->arg1;
    delete args;

    while (1){
        //std::pair<bam1_t *, bam1_t *> cur_read_pair;
        read_pair_t  cur_read_pair;
        while(1){
            int is_retrieved = 0;
            pthread_spin_lock(&thread_context->cur_reads_lock);
            
            if(thread_context->cur_reads->size()>0){
                cur_read_pair = thread_context->cur_reads->back();
                thread_context->cur_reads->pop_back();
                
                is_retrieved = 1;
            }
            pthread_spin_unlock(&thread_context->cur_reads_lock);
            if(global_context->all_finished && !is_retrieved) return NULL;
        
            if(is_retrieved) break;
            else usleep(tick_time);
        }
        process_aligned_fragment(global_context,thread_context,cur_read_pair);

    }
}


//start threads
void init_thread(global_context_t * global_context,unsigned short threadNumber,int mapq)
{
    global_context->thread_number = threadNumber;
    //global_context -> thread_contexts = (thread_context_t * ) malloc(sizeof(thread_context_t) * global_context -> thread_number);
   if(threadNumber >1){ 
    for(int i=0;i<threadNumber;i++)
    {
	    thread_context_t *th_contx = new thread_context_t();
        pthread_spin_init(&th_contx->cur_reads_lock, PTHREAD_PROCESS_PRIVATE);
		th_contx->thread_id = i;
       // th_contx->cur_reads = new std::vector<std::pair<bam1_t *,bam1_t *> >() ;
        
        th_contx->cur_reads = new std::vector<read_pair_t> ();
        
		th_contx->res = new Results();
        th_contx->clip_prof = new Clipping_prof(mapq);
        th_contx->cov_prof = new Coverage_prof(global_context->geneIdx);
        th_contx->inDist_prof = new InnerDist_prof();
        global_context->thread_contexts.push_back(th_contx);
    }
	
    for(int i=0;i<threadNumber;i++)
	{
		    struct arg_struct * args = new arg_struct();
			args->arg1 = global_context;
			args->arg2 = global_context->thread_contexts[i];
           //int ret = pthread_create(&global_context -> thread_contexts[i].thread_object, NULL, worker, (void *)&args);
           pthread_create(&global_context -> thread_contexts[i]->thread_object, NULL, worker, (void *)args);
   }
    }
	else {
	    thread_context_t *th_contx = new thread_context_t();
        //th_contx->cur_reads = new std::vector<std::pair<bam1_t *,bam1_t *> >() ;
        th_contx->cur_reads = new std::vector<read_pair_t> ();
        
		th_contx->res = new Results();
        th_contx->clip_prof = new Clipping_prof(mapq);
        th_contx->cov_prof = new Coverage_prof(global_context->geneIdx);
        th_contx->inDist_prof = new InnerDist_prof();
        global_context->thread_contexts.push_back(th_contx);
	  
	}

    return;
}

//release memory
void destroy_thread_context(global_context_t * global_context)
{
        
        for(int i=0; i < global_context-> thread_number; i++)
        {
            delete global_context->thread_contexts[i]->res;
            delete global_context->thread_contexts[i]->cov_prof;
            delete global_context->thread_contexts[i]->clip_prof;
            delete global_context->thread_contexts[i]->cur_reads;
            delete global_context->thread_contexts[i]->inDist_prof;
            if(global_context->thread_number >1){
            pthread_spin_destroy(&global_context -> thread_contexts[i]->cur_reads_lock);
			}
			delete global_context -> thread_contexts[i];
        }
    
}

//join threads
void join_threads(global_context_t * global_context)
{
    int xk1;
    for(xk1=0; xk1<global_context-> thread_number; xk1++)
        pthread_join(global_context -> thread_contexts[xk1]->thread_object, NULL);
}

//merge results
void merge_results(global_context_t * global_context,Results * res,Clipping_prof * clip_prof,Coverage_prof * cov_prof, InnerDist_prof * inDist_prof)
{
    for(int i=0;i<global_context->thread_number;i++)
    {
        res->add(global_context->thread_contexts[i]->res);
        clip_prof->add(global_context->thread_contexts[i]->clip_prof);
        cov_prof->add(global_context->thread_contexts[i]->cov_prof);
        inDist_prof->add(global_context->thread_contexts[i]->inDist_prof);
    }
}

// parsing BAM file
Results QC(std::string smp_name,GeneFeatures * geneIdx, rRNA * rRNAIdx,char * inputFile,std::string outdir,std::string outdir_fig,std::string stranded,int threadNumber,int mapq)
{
    //'''This class provides fuctions to parsing SAM or BAM files and quality controls.'''
    //'''constructor. input could be bam or sam'''

    int q_cut = mapq;
    global_context_t global_context;
    global_context.stranded = stranded;
    global_context.all_finished = 0;
    
    //open SAM/BAM file
    samFile *fp_in = NULL;
    bam1_t * aligned_read = NULL;
    bam_hdr_t *header = NULL;
    
    std::string format = "BAM";
    
    fp_in = sam_open(inputFile, "rb");
    
    if(NULL == fp_in) {
        std::cout << "Could not open file " << inputFile << std::endl;
        std::exit(1);
    }
    if(hts_get_format(fp_in)->format == sam) {format = "SAM";}
    global_context.format = format;
    
    aligned_read = bam_init1();
    
    Results main_res;
    //read header
    header = sam_hdr_read(fp_in);
	if(NULL == header){
        std::cout << "No header information." << std::endl;
		std::exit(1);
    }
    for(int i=0;i< header->n_targets;i++)
	{
	   char * p = header->target_name[i];
	   std::string chr(p);
	   global_context.refnames.push_back(chr);
	}
	//return main_res;

    std::string prev_read_id = "";
    std::vector<bam1_t *> multi_read1;
    std::vector<bam1_t *> multi_read2;
    
    std::vector<bam1_t *> alignments_per_read;
    

    Clipping_prof main_clip_prof(outdir+smp_name,outdir_fig+smp_name,mapq);
    Coverage_prof main_cov_prof(outdir+smp_name,outdir_fig+smp_name,geneIdx);
    //ReadDup_prof main_rDup_prof(outdir,outdir_fig);
    InnerDist_prof main_inDist_prof(outdir+smp_name,outdir_fig+smp_name);

    int lineno = 0;
	global_context.geneIdx = geneIdx;
	global_context.rRNAIdx = rRNAIdx;
    init_thread(&global_context, threadNumber,mapq);
   
    int current_thread_id = 0;
    
    while(sam_read1(fp_in, header,aligned_read) > 0)
    {
        lineno +=1 ;
        //std::cout << lineno << std::endl;
        
        std::string qname ="";
        char * name = bam_get_qname(aligned_read);
        if(name !=NULL)
            {
                qname = std::string(name);
            }
            //SE reads
            if (! IS_PAIRED(aligned_read))
            {
                //count unmapped read
                if (IS_UNMAPPED(aligned_read))
                {
                    main_res.unmapped_reads += 1;
                    main_res.total_reads += 1;
                    bam_destroy1(aligned_read);
                    aligned_read =  bam_init1();
                    continue;
                }
                
                if (qname == prev_read_id || prev_read_id == "") {
                    alignments_per_read.push_back(aligned_read);
                    prev_read_id = qname;
                    aligned_read =  bam_init1();
                    continue;
                }
                else {
                    bam1_t *cur_read = NULL;
                    if (alignments_per_read.size() == 1){  //unique read
                        main_res.uniq_mapped_reads += 1;
                        main_res.total_reads += 1;
                        
                        cur_read = alignments_per_read[0];
                        
                        int skip_read  = 0;
                        
                        if (cur_read->core.flag & BAM_FQCFAIL || cur_read->core.qual < q_cut){    //skip QC fail read
                            main_clip_prof.set_qual(cur_read->core.qual);
                            main_res.low_qual += 1;
                            skip_read = 1;
                        }
                        if (IS_DUP(cur_read)){       //skip duplicate read
                            skip_read = 1;
                        }
                        if (skip_read == 1) {
                            bam_destroy1(cur_read);
                            cur_read = NULL;
                        }
                        if (cur_read ){
                            if (IS_REVERSE(cur_read)) {main_res.reverse_read += 1;}
                            else {main_res.forward_read += 1;}
                            read_pair_t read_pair; // = new read_pair_t();
                            read_pair.first = cur_read;
                            read_pair.second = NULL;
                            read_pair.type = 0;
                            
                            if(global_context.thread_number >1){

                                thread_context_t * thread_context = global_context.thread_contexts[current_thread_id];
                                pthread_spin_lock(&thread_context->cur_reads_lock);
                                
                                //thread_context->cur_reads->push_back(std::pair<bam1_t *,bam1_t *>(cur_read,NULL));
                                thread_context->cur_reads->push_back(read_pair);
                                
                                pthread_spin_unlock(&thread_context->cur_reads_lock);
                                current_thread_id++;
                                if(current_thread_id >= global_context.thread_number) current_thread_id = 0;
                            }
                            else {
                               // process_aligned_fragment(&global_context,global_context.thread_contexts[0],std::pair<bam1_t *,bam1_t *> (cur_read,NULL));
                                process_aligned_fragment(&global_context,global_context.thread_contexts[0],read_pair);
                            } 
                        }

                    }
                    else { //multi reads
                        main_res.multi_mapped_reads += 1;
                        main_res.total_reads += 1;
                        cur_read = alignments_per_read[0];
                        //check whether it's rRNA read using only one of the multialignments
                        main_clip_prof.set_qual(cur_read->core.qual);
                        
                        //for(unsigned int i=0;i<alignments_per_read.size();i++){
                        for(unsigned int i=1;i<alignments_per_read.size();i++){
                            bam_destroy1(alignments_per_read[i]);
                        }
                        //cur_read = NULL;
                        read_pair_t  read_pair; // = new read_pair_t();
                        read_pair.first = cur_read;
                        read_pair.second = NULL;
                        read_pair.type = 1;
                        
                        if(global_context.thread_number >1){
                        //thread_context_t * thread_context = global_context.thread_contexts + current_thread_id;
                        thread_context_t * thread_context = global_context.thread_contexts[current_thread_id];
                        pthread_spin_lock(&thread_context->cur_reads_lock);
                            
                        thread_context->cur_reads->push_back(read_pair);
                        //thread_context->cur_reads->push_back(std::pair<bam1_t *,bam1_t *>(cur_read,NULL));
                        pthread_spin_unlock(&thread_context->cur_reads_lock);
                        current_thread_id++;
                        if(current_thread_id >= global_context.thread_number) current_thread_id = 0;
                        }
						else {
                           //process_aligned_fragment(&global_context,global_context.thread_contexts[0],std::pair<bam1_t *,bam1_t *> (cur_read,NULL));
                            process_aligned_fragment(&global_context,global_context.thread_contexts[0],read_pair);
                        
                       } 
                    }
                    alignments_per_read.clear();
                    alignments_per_read.push_back(aligned_read);
                    prev_read_id = qname;
					aligned_read = bam_init1();
                }
            }
            else {//pair end read
                main_res.is_pairEnd = true;
                size_t flag_pos = qname.find('/');
                //cur_read_id = aligned_read.qname;
                if (flag_pos != std::string::npos) {qname = qname.substr(0,flag_pos);}

                if (qname == prev_read_id ){
                    if (IS_READ1(aligned_read) ){ multi_read1.push_back(aligned_read);}
                    if (IS_READ2(aligned_read) ){multi_read2.push_back(aligned_read);}
					aligned_read = bam_init1();
					continue;
                }
                else {
                    bam1_t * cur_read1 = NULL;
                    bam1_t * cur_read2 = NULL;
                    main_res.total_reads += 1;
                        
                    //multi-reads
                    if (multi_read1.size() >1 || multi_read2.size() > 1){
                        main_res.multi_mapped_reads += 1;
                        read_pair_t  read_pair;
                        read_pair.first = NULL;
                        read_pair.second = NULL;
                        read_pair.type = 1;
                        
                        if (multi_read1.size() > 1 ){
                            main_clip_prof.set_qual(multi_read1[0]->core.qual);
                            read_pair.first = multi_read1[0];
                        }
                        if (multi_read2.size() > 1 ){
                            main_clip_prof.set_qual(multi_read2[0]->core.qual);
                            read_pair.second = multi_read2[0];
                        }
                        //for(unsigned int i=0;i< multi_read1.size();i++)
                        for(unsigned int i=1;i< multi_read1.size();i++)
                        {
                            bam_destroy1(multi_read1[i]);
                        }
                        //for(unsigned int i=0;i< multi_read2.size();i++)
                        for(unsigned int i=1;i< multi_read2.size();i++)
                        {
                            bam_destroy1(multi_read2[i]);
                        }
                        //read_pair_t * read_pair = new read_pair_t();
                        
                        
                        if(global_context.thread_number > 1){
                            //thread_context_t * thread_context = global_context.thread_contexts + current_thread_id;
                            thread_context_t * thread_context = global_context.thread_contexts[current_thread_id];
                            
                            pthread_spin_lock(&thread_context->cur_reads_lock);
                            
                            //thread_context->cur_reads->push_back(std::pair<bam1_t *,bam1_t *>(cur_read1,cur_read2));
                            thread_context->cur_reads->push_back(read_pair);
                            
                            pthread_spin_unlock(&thread_context->cur_reads_lock);
                            current_thread_id++;
                            if(current_thread_id >= global_context.thread_number)
                                current_thread_id = 0;
                        }
                        else {
                            //process_aligned_fragment(&global_context,global_context.thread_contexts[0],std::pair<bam1_t *,bam1_t *> (cur_read1,cur_read2));
                            process_aligned_fragment(&global_context,global_context.thread_contexts[0],read_pair);
                        }
                    }
                    else {
                        //uniq read
                        if (multi_read1.size() == 1) {cur_read1 = multi_read1[0];}
                        if (multi_read2.size() == 1){ cur_read2 = multi_read2[0];}

                        if (cur_read1 && !IS_UNMAPPED(cur_read1)){
                            
                            if ( cur_read1->core.flag & BAM_FQCFAIL || cur_read1->core.qual < q_cut ){
                                main_clip_prof.set_qual(cur_read1->core.qual);
                                bam_destroy1(cur_read1);
                                cur_read1 = NULL;
                                main_res.low_qual_read1 += 1;
                            }
                            else {
                            main_res.mapped_read1 += 1;
                                if (IS_REVERSE(cur_read1)){ main_res.reverse_read += 1;}
                                else{ main_res.forward_read +=1;}
                                }
                        }
                        else {
                            if(cur_read1) { bam_destroy1(cur_read1); cur_read1 = NULL;}
                            if (prev_read_id != ""){ main_res.unmapped_read1 += 1; }
                        }
                            
                        if (cur_read2 && !IS_UNMAPPED(cur_read2)){
                            
                            if (cur_read2->core.flag & BAM_FQCFAIL || cur_read2->core.qual < q_cut){
                                main_clip_prof.set_qual(cur_read2->core.qual);
                                main_res.low_qual_read2 += 1;
                                bam_destroy1(cur_read2);
                                cur_read2 = NULL;
                            }
                            else {
                            main_res.mapped_read2 += 1;
                                if (!cur_read1){
                                    if (!IS_REVERSE(cur_read2)){ main_res.reverse_read += 1;}
                                    else{ main_res.forward_read +=1;}
                                }
                            }
                        }
                        else {
                            if(cur_read2) {bam_destroy1(cur_read2); cur_read2 = NULL;}
                            if (prev_read_id != ""){ // #not the first fragment
                                main_res.unmapped_read2 += 1;}
                        }
                            
                        if (cur_read1  || cur_read2 ){
                            //read_pair_t * read_pair = new read_pair_t();
                            read_pair_t  read_pair;
                            read_pair.first = cur_read1;
                            read_pair.second = cur_read2;
                            read_pair.type = 0;
                            
						if(global_context.thread_number > 1){
                            //thread_context_t * thread_context = global_context.thread_contexts + current_thread_id;
                            thread_context_t * thread_context = global_context.thread_contexts[current_thread_id];
                            
							pthread_spin_lock(&thread_context->cur_reads_lock);
                            
                            //thread_context->cur_reads->push_back(std::pair<bam1_t *,bam1_t *>(cur_read1,cur_read2));
                            thread_context->cur_reads->push_back(read_pair);
                            
                            pthread_spin_unlock(&thread_context->cur_reads_lock);
                            current_thread_id++;
                            if(current_thread_id >= global_context.thread_number) 
							    current_thread_id = 0;
                            }
							else {
                              //process_aligned_fragment(&global_context,global_context.thread_contexts[0],std::pair<bam1_t *,bam1_t *> (cur_read1,cur_read2));
                              process_aligned_fragment(&global_context,global_context.thread_contexts[0],read_pair);
                            }
                        }
                    }
                    multi_read1.clear();
                    multi_read2.clear();
                    prev_read_id = qname;
                    if (IS_READ1(aligned_read)){
                        multi_read1.push_back(aligned_read);}
                    if (IS_READ2(aligned_read)){
                        multi_read2.push_back(aligned_read);}
                    }
                
                aligned_read =  bam_init1();
            }
    }
    bam_destroy1(aligned_read);
	bam_hdr_destroy(header);
    sam_close(fp_in);
    global_context.all_finished = 1;
    //join threads
    if(global_context.thread_number >1)
    {
        join_threads(&global_context);
		}
        merge_results(&global_context,&main_res,&main_clip_prof,&main_cov_prof,&main_inDist_prof);
        destroy_thread_context(&global_context);
    
    
    main_res.splice = int(main_res.splice);
    main_res.noSplice = int(main_res.noSplice);
     
    if (main_res.is_pairEnd ){
        main_clip_prof.write(main_res.total_reads*2-main_res.unmapped_read1-main_res.unmapped_read2);
    }
    else {
        main_clip_prof.write(main_res.total_reads);
    }

    //main_rDup_prof.write();
    main_inDist_prof.write();
    int zero_exons = main_cov_prof.write(main_res.total_reads);

    main_res.read_dist_plot_file1 = outdir_fig +smp_name+ ".read_distr.png";
    main_res.read_dist_plot_file2 = outdir_fig + smp_name + ".read_distr_pie.png";
     
    //std::cout << main_res.rRNA_read << std::endl;   
    try {
            std::ofstream ROUT;
        
            ROUT.open (outdir +smp_name+ ".read_distr.r", std::ofstream::out);
            //ROUT = fopen(outfile + '.read_distr.r', 'w')
            ROUT << "png(\"" << main_res.read_dist_plot_file1 << "\",width=500,height=500,units=\"px\")\n";
        ROUT << "M=c(" << std::to_string(main_res.cds_exon_read) << "," << std::to_string(main_res.utr_5_read) << "," << std::to_string(main_res.utr_3_read) << "," << std::to_string(main_res.intron_read) << "," << std::to_string(main_res.intergenic_up1kb_read) << "," << std::to_string(main_res.intergenic_down1kb_read) << "," << std::to_string(main_res.rRNA_read) << "," << std::to_string(main_res.intergenic_read) << ")\n";
            
            ROUT << "Mname=c(\"CDS\",\"5UTR\",\"3UTR\",\"Intron\",\"TSS_Up_1Kb\",\"TES_Down_1Kb\",\"rRNA\",\"Others\")\n";
            ROUT << "val = barplot(M,xlab=\"\",space=1,ylab=\"Read Counts\",col=\"blue\",border=\"NA\")\n";
            ROUT << "text(x=seq(val[1],val[8],by=2),y=rep(0,8),srt=60,adj=0,offset=2,pos=1,xpd=T,labels=Mname)\n";
            ROUT << "dev.state = dev.off()\n";
        ROUT.close();

            ROUT.open(outdir + smp_name+".read_distr_pie.r", std::ofstream::out);
            if (geneIdx->total_exon != 0 ){
                ROUT << "png(\"" << main_res.read_dist_plot_file2 << "\",width=500,height=500,units=\"px\")\n";
                ROUT << "pie(c(" << std::to_string(geneIdx->total_exon-zero_exons) << ',' << std::to_string(zero_exons) << "),labels=c(\"Covered  " << (geneIdx->total_exon - zero_exons) << " exons\",\"Uncovered\"),main=\"Exons\",radius=0.6,clockwise=T,col=c(\"blue\",\"white\"))\n";
                ROUT << "dev.state = dev.off()\n";
            }
        ROUT.close();
    }catch(std::ofstream::failure e ){
            std::cout << "Error in writing plotting scripts.\n" << std::endl;
    }

    main_res.insert_plot_file = main_inDist_prof.InnDist_fig_file;
    main_res.insert_file = main_inDist_prof.InnDist_data_file;
    
    main_res.clipping_plot_file = main_clip_prof.clip_fig_file;
    main_res.mapq_plot_file = main_clip_prof.mapq_fig_file;
    main_res.mapq_file = main_clip_prof.mapq_data_file;
    
    //res.read_dup_plot_file = rDup_prof.plot_file
    main_res.readLen_plot_file = main_clip_prof.readlen_fig_file;
    
    main_res.read_cov_plot_file = main_cov_prof.cov_fig_file;
    main_res.geneCount_file = main_cov_prof.transcov_data_file;
    main_res.trans_cov_plot_file = main_cov_prof.transcov_fig_file;
        //res.seqDeDup_percent = rDup_prof.seqDeDup_percent
        //res.posDeDup_percent = rDup_prof.posDeDup_percent
    return main_res;
}

int run_qc(char* out_dir, char* outfig_dir,char * ann_file, char* attrID, char* input_file,char* rRNA_file,char* label,int mapq,char* stranded,int thread_num)
{
    std::string gtf_fname (ann_file);
    std::string id_attrID (attrID);
    //std::string ifile (input_file);
    std::string smp_name (label);
    std::string strand_info (stranded);
    std::string data_outdir (out_dir);
    std::string fig_outdir (outfig_dir);
    std::string smp_res_fname = data_outdir+smp_name+".res.txt";
    
    //int thread_num = 1;
    
    std::string rRNA_fname (rRNA_file);
    rRNA * rRNAIdx = NULL;
    
    GeneFeatures * geneIdx = new GeneFeatures(gtf_fname,id_attrID);
    if (rRNA_fname != "") {
        rRNAIdx = new rRNA(rRNA_fname);
    }
    Results res = QC(smp_name,geneIdx,rRNAIdx,input_file,data_outdir,fig_outdir,strand_info,thread_num,mapq);
    
    res.write(smp_res_fname);
    
    delete geneIdx;
    delete rRNAIdx;
    
    return 1;
}



