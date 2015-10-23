#!/usr/bin/env python
'''data structure store results.'''



class Results :
    
    def __init__(self):
        self.filename = ""
        self.is_pairEnd = False
        self.clipping_plot_file = ""
        self.mapq_plot_file = ""
        self.mapq_file = ""
        self.read_cov_plot_file = ""
        self.trans_cov_plot_file = ""
        self.insert_plot_file = ""
        self.insert_file = ""
        self.read_dist_plot_file1 = ""
        self.read_dist_plot_file2 = ""
        self.read_dup_plot_file = ""
        self.readLen_plot_file = ""
        self.geneCount_file = ""
        
        self.seqDeDup_percent = 0
        self.posDeDup_percent = 0
        
        self.no_clipping = False
        self.no_rRNA = False
        
        
        self.total_reads = 0
        self.uniq_mapped_reads = 0
        self.multi_mapped_reads = 0
        self.unmapped_reads = 0
        self.low_qual = 0
        self.low_qual_read1 = 0
        self.low_qual_read2 = 0
        self.pcr_dup = 0
        
        self.unmapped_read1 = 0
        self.unmapped_read2 = 0
        self.mapped_read1 = 0
        self.mapped_read2 = 0 
        self.forward_read = 0
        self.reverse_read = 0
        self.paired_reads = 0
        
        self.mapped_plus_minus = 0
        self.mapped_plus_plus = 0
        self.mapped_minus_plus = 0
        self.mapped_minus_minus = 0
        
        self.ins_read = 0
        self.del_read = 0
        
        self.noSplice = 0
        self.splice = 0
        self.paired_diff_chrom = 0


