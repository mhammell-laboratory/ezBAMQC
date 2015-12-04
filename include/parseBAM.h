//
//  parseBAM.h
//  BAMQC_c++
//
//  Created by Ying Jin on 11/18/15.
//  Copyright (c) 2015 Ying Jin. All rights reserved.
//

#ifndef __BAMQC_c____parseBAM__
#define __BAMQC_c____parseBAM__

#include <stdio.h>

extern "C" {
int run_qc(char* out_dir, char* outfig_dir,char * ann_file, char* attrID, char* input_file,char* rRNA_file,char* label,int maqp,char* stranded,int threadNum);
}
#endif /* defined(__BAMQC_c____parseBAM__) */
