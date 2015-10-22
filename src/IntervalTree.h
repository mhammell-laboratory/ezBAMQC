//
//  IntervalTree.h
//  BAMQC-0.5
//
//  Created by Ying Jin on 9/15/15.
//  Copyright (c) 2015 Ying Jin. All rights reserved.
//

#ifndef __BAMQC_0_5__IntervalTree__
#define __BAMQC_0_5__IntervalTree__

#include <stdio.h>
#include <string>
#include <vector>



class Interval{
public:
    long  start;
    long stop;
    std::string gene;
    std::string type;

    Interval();
    Interval(std::string g, long st, long end, std::string t);
    ~Interval();
};





class IntervalTree{
public :
    std::vector<Interval> itvlist;
    IntervalTree* left ;
    IntervalTree* right ;
    float center ;

    IntervalTree();
    IntervalTree(std::vector<Interval> intervals, int depth=16, int minbucket=16, long extent_st=-1, long extent_end = -1 , int maxbucket=512);
    ~IntervalTree();
    std::vector<Interval> find(long start, long stop);
    std::vector<std::string>  find_gene(long start, long stop);
private:
    void clear(IntervalTree * node);

};



#endif /* defined(__BAMQC_0_5__IntervalTree__) */
