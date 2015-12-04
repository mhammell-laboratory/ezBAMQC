//
//  Constants.h
//  BAMQC_c++
//
//  Created by Ying Jin on 10/28/15.
//  Copyright (c) 2015 Ying Jin. All rights reserved.
//

#ifndef BAMQC_c___Constants_h
#define BAMQC_c___Constants_h

#include <limits>

#define BIN_SIZE  100000

#define MAX_BUCKET  128
#define MIN_BUCKET  16
#define DEPTH  16

//#define SAMPLESIZE 500000
#define SAMPLESIZE std::numeric_limits<int>::max()
#define LOW_BOUND -200
#define UPPER_BOUND 1000
#define STEP 10

#define CDS 1
#define UTR5 2
#define UTR3 3
#define INTRON 4
#define ITGUP1K 5
#define ITGDN1K 6
#define INTERGENIC 7
#define RRNA 8

#define MAX_READ_LEN 1000

#endif
