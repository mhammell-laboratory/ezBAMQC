"""Module Description

    Copyright (c) 2014, Ying Jin <yjin@cshl.edu >


    This code is free software; you can redistribute it and/or modify it
    under the terms of the Artistic License (see the file COPYING included
    with the distribution).

    @author:  Ying Jin
    @contact: yjin@cshl.edu
    """
import sys, time,re
import logging
import gzip
from math import ceil,floor
import collections

from libBAMQC.IntervalTree import *


class GFF_Reader( ):

   """Parse a GFF file

   Pass the constructor either a file name or an iterator of lines of a
   GFF files. If a file name is specified, it may refer to a gzip compressed
   file.

   Yields tuple of (gene_id,chrom,strand,start position,end position,type)

   """

   def __init__( self, filename, id_attribute):
      self.line_no = None
      self.filename = filename
      self.id_attribute = id_attribute
      self._re_attr_main = re.compile( "\s*([^\s\=]+)[\s=]+(.*)" )
   

   def __iter__( self ):
      self.line_no = 0
      if self.filename.lower().endswith( ( ".gz" , ".gzip" ) ):
            lines = gzip.open( self.filename )
      else:
            lines = open( self.filename )

      for line in lines:
          self.line_no += 1
          if line == "\n" or line.startswith('#'):
              continue
          ( seqname, source, feature, start, end, score, strand, frame, attributeStr ) = line.split("\t")
          id = self.__parse_GFF_attr_string(attributeStr,self.id_attribute)

          yield (id, seqname, strand, int(start), int(end), feature)

      lines.close()
      self.line_no = None

   def __parse_GFF_attr_string(self,attributeStr,id_interested) :


       for pairs in attributeStr.split(';') :
           if pairs.count('"') not in [0,2] :
               raise ValueError, "The attribute string seems to contain mismatched quotes."
           nv = self._re_attr_main.match(pairs)
           if not nv :
               raise ValueError, "Failure parsing GFF attribute line."
           val = nv.group(2)
           name = nv.group(1)
           if name == id_interested :
               return val
       return None

   def get_line_number_string( self ):
      if self.line_no is None:
            return "file %s closed" % self.filename

      else:
         return "line %d of file %s" % ( self.line_no, self.filename )


class GeneFeatures:
    """index of Gene annotations.
        """
    def __init__ (self):
        #self.featureIdxs_plus = {}
        #self.featureIdxs_minus = {}
        #self.featureIdxs_nostrand = {}
        self.features = []
        #gene as key , exons as values
        self.plus_gene_exons = {}
        self.minus_gene_exons = {}
        
        self.cds_exon_idx_plus = {}
        self.cds_exon_idx_minus = {}
        #self.cds_exon_idx_nostrand = {}
        
        self.intron_idx_plus = {}
        self.intron_idx_minus = {}
        #self.intron_idx_nostrand = {}
        
        self.utr_3_idx_plus  = {}
        self.utr_3_idx_minus = {}
        #self.utr_3_idx_nostrand = {}
        
        self.utr_5_idx_plus = {}
        self.utr_5_idx_minus = {}
        #self.utr_5_idx_notrand = {}
        
        self.intergenic_up_1kb_plus = {}
        self.intergenic_up_1kb_minus = {}
        #self.intergenic_up_1kb_nostrand = {}
        
        self.intergenic_down_1kb_plus = {}
        self.intergenic_down_1kb_minus = {}
        #self.intergenic_down_1kb_nostrand = {}
        
    def set(self,GTFfilename,id_attribute) :
        self.read_features(GTFfilename,id_attribute)


    # Reading & processing annotation files
    def read_features(self,gff_filename, id_attribute) :

        #dict of dicts since the builtin type doesn't support it for some reason
        temp_plus = collections.defaultdict(dict)
        temp_minus = collections.defaultdict(dict)
        #temp_nostrand = collections.defaultdict(dict)
        
        cds_plus = dict()
        cds_minus = dict()
        #cds_nostrand = dict()
            
        intron_plus = dict()
        intron_minus = dict()
        #intron_nostrand = dict()
            
        utr5_plus = dict()
        utr5_minus = dict()
        #utr5_nostrand = dict()
            
        utr3_plus = dict()
        utr3_minus = dict()
        #utr3_nostrand = dict()
        
        itgUp1k_plus = dict()
        itgUp1k_minus = dict()
        #itgUp1k_nostrand = dict()
        
        itgDn1k_plus = dict()
        itgDn1k_minus = dict()
        #itgDn1k_nostrand = dict()

        # read count of features in GTF file
        gff = GFF_Reader(gff_filename,id_attribute)  # (id, seqname, strand, int(start), int(end), feature)
        i = 0
        counts = 0
        try:
            for (id, chrom, strand, start, end, feature) in gff:
                if id is None :
                    continue
                #save gene id
                if id not in self.features :
                    self.features.append(id)
                
                st_end_str = str(start)+":"+str(end)
                if feature == "CDS" :
                    
                    if strand == "+" :
                        if id in self.plus_gene_exons :
                            if st_end_str not in self.plus_gene_exons[id] :
                                self.plus_gene_exons[id].append(st_end_str)
                                
                                try:
                                    cds_plus[chrom].append(Interval(id,start,end))
                                except:
                                    cds_plus[chrom] = [Interval(id,start,end)]
                        else :
                            self.plus_gene_exons[id] = [st_end_str]
                            try:
                                cds_plus[chrom].append(Interval(id,start,end))
                            except:
                                cds_plus[chrom] = [Interval(id,start,end)]
                
                    if strand == "-" :
                        if id in self.minus_gene_exons :
                            if st_end_str not in self.minus_gene_exons[id] :
                                self.minus_gene_exons[id].append(st_end_str)
                    
                        else :
                            self.minus_gene_exons[id] = [st_end_str]
                            try:
                                cds_minus[chrom].append(Interval(id,start,end))
                            except:
                                cds_minus[chrom] = [Interval(id,start,end)]
                                        
                                        
                if feature == "exon":
                    counts += 1
                    
                    if strand == "+"  :
                        if id in self.plus_gene_exons :
                            if st_end_str not in self.plus_gene_exons[id] :
                                self.plus_gene_exons[id].append(st_end_str)
                                try:
                                    temp_plus[chrom][id].append((start,end))
                                except:
                                    temp_plus[chrom][id] = [(start,end)]
                        else :
                            self.plus_gene_exons[id] = [st_end_str]
                            try:
                                temp_plus[chrom][id].append((start,end))
                            except:
                                temp_plus[chrom][id] = [(start,end)]
    
                    if strand == "-" :
                        if id in self.minus_gene_exons :
                            if st_end_str not in self.minus_gene_exons[id] :
                                self.minus_gene_exons[id].append(st_end_str)
                                try:
                                    temp_minus[chrom][id].append((start,end))
                                except:
                                    temp_minus[chrom][id] = [(start,end)]
                        else :
                            self.minus_gene_exons[id] = [st_end_str]
                            try:
                                temp_minus[chrom][id].append((start,end))
                            except :
                                temp_minus[chrom][id] = [(start,end)]
                                    
                i += 1
                if i % 100000 == 0 :
                    sys.stderr.write("%d GTF lines processed.\n" % i)
        except:
            sys.stderr.write("Error occured in %s.\n" % gff.get_line_number_string())
            raise

        if counts == 0 :
            sys.stderr.write("Warning: No features of type 'exon' or 'CDS' found in gene GTF file.\n")

        #build interval trees
        #trees for CDS
        #sys.stderr.write("start to build interval tree\n")
        for chr in cds_plus.keys() :
            self.cds_exon_idx_plus[chr] = IntervalTree(cds_plus[chr])
        
        #sys.stderr.write("done with cds exon plus interval tree\n")
        for chr in cds_minus.keys() :
            self.cds_exon_idx_minus[chr] = IntervalTree(cds_minus[chr])
        
        
        #sys.stderr.write("start to build exon interval tree\n")
        #trees for the others
        for each_chrom in temp_plus.keys():
            inputlist = []
            intronlist = []
            utr5list = []
            utr3list = []
            itgUp1klist = []
            itgDn1klist = []
          
            for each_gene in temp_plus[each_chrom].keys():
                #sort exons
                #sys.stderr.write(each_chrom + "\t" + each_gene + "\n")
                sorted_exons = sorted(temp_plus[each_chrom][each_gene],key=lambda tup:tup[0])
                #sys.stderr.write("sorted_exon\n")
                #exon start/end position
                exon_starts =  map(lambda tup: tup[0],sorted_exons)
                exon_ends =  map(lambda tup: tup[1],sorted_exons)
                
                #sys.stderr.write("sorted_exon_ends\n")
                #intron start/end position
                intron_starts = map(lambda x: x + 1, exon_ends[:-1])
                intron_ends =map(lambda x: x - 1, exon_starts[1:])
                
                #sys.stderr.write("sorted_intron_ends\n")
                #sys.stderr.write("exon starts length " + str(len(exon_starts)) + "\n")
                #5/3 UTR position
                if len(exon_starts) > 1 :
                    utr5_st = exon_starts[0]
                    utr5_end = min(exon_ends[0],exon_starts[1])
                    try :
                        utr5list.append(Interval(each_gene,utr5_st,utr5_end))
                    except:
                        utr5list = [Interval(each_gene,utr5_st,utr5_end)]
                    
                    utr3_st = max(exon_starts[-1],exon_ends[-2])
                    utr3_end = exon_ends[-1]
                    try:
                        utr3list.append(Interval(each_gene,utr3_st,utr3_end))
                    except:
                        utr3list = [Interval(each_gene,utr3_st,utr3_end)]
                itgUp1k_st = max(0,exon_starts[0]-1000)
                itgUp1k_end = max(0,exon_starts[0]-1 )
                
                try:
                    itgUp1klist.append(Interval(each_gene,itgUp1k_st,itgUp1k_end))
                except:
                    itgUp1klist = [Interval(each_gene,itgUp1k_st,itgUp1k_end)]
                
                itgDn1k_st = exon_ends[-1]+1
                itgDn1k_end = exon_ends[-1] + 1000
                
                try:
                    itgDn1klist.append(Interval(each_gene,itgDn1k_st,itgDn1k_end))
                except:
                    itgDn1klist = [Interval(each_gene,itgDn1k_st,itgDn1k_end)]
                
                #sys.stderr.write("before_intron_lists\n")
                
                for st,end in zip(intron_starts,intron_ends):
                    if st < end - 1 :
                        try:
                            intronlist.append(Interval(each_gene,st,end))
                        except:
                            intronlist = [Interval(each_gene,st,end)]
                if len(cds_plus) + len(cds_minus) == 0 :
                    for (start,end) in temp_plus[each_chrom][each_gene]:
                        try:
                            inputlist.append(Interval(each_gene,start,end))
                        except:
                            inputlist = [Interval(each_gene,start,end)]
            if len(inputlist) > 0 :
                self.cds_exon_idx_plus[each_chrom] = IntervalTree(inputlist)
            if len(intronlist) > 0 :
                self.intron_idx_plus[each_chrom] = IntervalTree(intronlist)
            if len(utr3list) > 0 :
                self.utr_3_idx_plus[each_chrom] = IntervalTree(utr3list)
            if len(utr5list) > 0 :
                self.utr_5_idx_plus[each_chrom] = IntervalTree(utr5list)
            if len(itgUp1klist) > 0 :
                self.intergenic_up_1kb_plus[each_chrom] = IntervalTree(itgUp1klist)
            if len(itgDn1klist) > 0 :
                self.intergenic_down_1kb_plus[each_chrom] = IntervalTree(itgDn1klist)
                            
                            #sys.stderr.write("done with exon plus interval tree\n")
        #trees for exons
        for each_chrom in temp_minus.keys():
            inputlist = []
            intronlist = []
            utr5list = []
            utr3list = []
            itgUp1klist = []
            itgDn1klist = []

            for each_gene in temp_minus[each_chrom].keys():
                #sort exons
                sorted_exons = sorted(temp_minus[each_chrom][each_gene],key=lambda tup:tup[0])
                #exon start/end position
                exon_starts =  map(lambda tup: tup[0],sorted_exons)
                exon_ends =  map(lambda tup: tup[1],sorted_exons)
                
                #intron start/end position
                intron_starts = map(lambda x: x + 1, exon_ends[:-1])
                intron_ends =map(lambda x: x - 1, exon_starts[1:])
                
                #5/3 UTR position
                if len(exon_starts) > 1 :
                    utr3_st = exon_starts[0]
                    utr3_end = min(exon_ends[0],exon_starts[1])
                    
                    try:
                        utr3list.append(Interval(each_gene,utr3_st,utr3_end))
                    except:
                        utr3list = [Interval(each_gene,utr3_st,utr3_end)]
                    
                    utr5_st = max(exon_starts[-1],exon_ends[-2])
                    utr5_end = exon_ends[-1]
                    try:
                        utr5list.append(Interval(each_gene,utr5_st,utr5_end))
                    except:
                        utr5list = [Interval(each_gene,utr5_st,utr5_end)]
                        
                
                itgDn1k_st = max(0,exon_starts[0]-1000)
                itgDn1k_end = max(0,exon_starts[0]-1)
                
                try:
                    itgUp1klist.append(Interval(each_gene,itgUp1k_st,itgUp1k_end))
                except:
                    itgUp1klist = [Interval(each_gene,itgUp1k_st,itgUp1k_end)]
                
                itgUp1k_st = exon_ends[-1]+1
                itgUp1k_end = exon_ends[-1] + 1000
                
                try:
                    itgDn1klist.append(Interval(each_gene,itgDn1k_st,itgDn1k_end))
                except:
                    itgDn1klist = [Interval(each_gene,itgDn1k_st,itgDn1k_end)]
                
                for st,end in zip(intron_starts,intron_ends):
                    if st < end - 1 :
                        try:
                            intronlist.append(Interval(each_gene,st,end))
                        except:
                            intronlist = [Interval(each_gene,st,end)]
                
                if len(cds_plus) + len(cds_minus) == 0 :
                    for (start,end) in temp_minus[each_chrom][each_gene]:
                        try:
                            inputlist.append(Interval(each_gene,start,end))
                        except:
                            inputlist = [Interval(each_gene,start,end)]
            if len(inputlist) > 0 :
                self.cds_exon_idx_minus[each_chrom] = IntervalTree(inputlist)
            
            if len(intronlist) > 0 :
                self.intron_idx_minus[each_chrom] = IntervalTree(intronlist)
            if len(utr3list) > 0 :
                self.utr_3_idx_minus[each_chrom] = IntervalTree(utr3list)
            if len(utr5list) > 0 :
                self.utr_5_idx_plus[each_chrom] = IntervalTree(utr5list)
            if len(itgUp1klist) > 0 :
                self.intergenic_up_1kb_minus[each_chrom] = IntervalTree(itgUp1klist)
            if len(itgDn1klist) > 0 :
                self.intergenic_down_1kb_minus[each_chrom] = IntervalTree(itgDn1klist)


    #find exons of given gene that overlap with the given intervals
    #return list of tuples 
    def find_exons(self,chrom,st,end,strand,gene) :
        exons = []
        fs = []
        try:
            if strand == "+" or strand == ".":
                if chrom in self.cds_exon_idx_plus :
                    fs = self.cds_exon_idx_plus[chrom].find(st,end)
                
                if chrom in self.utr_5_idx_plus:
                    fs += self.utr_5_idx_plus[chrom].find(st,end)
                if chrom in self.utr_3_idx_plus :
                    fs += self.utr_3_idx_plus[chrom].find(st,end)
            
            if strand == "-" or strand == ".":
                if chrom in self.cds_exon_idx_minus:
                    fs += self.cds_exon_idx_minus[chrom].find(st, end)

                if chrom in self.utr_5_idx_minus:
                    fs = self.utr_5_idx_minus[chrom].find(st, end)
                if chrom in self.utr_3_idx_minus :
                    fs += self.utr_3_idx_minus[chrom].find(st,end)
                   
        except:
            raise
        
        for itv in fs :
            if itv.gene == gene :
                exons.append((itv.start,itv.stop))
        
        exons = sorted(exons,key=lambda tup:tup[0])
        return exons
    
    def ovp_utr5(self,itv_list,strand) :
        fs = []
        genes = []
        for itv in itv_list :
            try:
                if strand == "+" or strand == ".":
                    if itv[0] in self.utr_5_idx_plus:
                        fs = self.utr_5_idx_plus[itv[0]].find_gene(itv[1],itv[2])
        
                if strand == "-" or strand == ".":
                    if itv[0] in self.utr_5_idx_minus:
                        fs = fs + self.utr_5_idx_minus[itv[0]].find_gene(itv[1], itv[2])
                if len(fs) > 0:
                    genes = genes + fs
        
            except:
                raise
        genes = sorted(set(genes))
        return genes
    #
    def ovp_utr3(self,itv_list,strand) :
        fs = []
        genes = []
        for itv in itv_list :
            try:
                if strand == "+" or strand == ".":
                    if itv[0] in self.utr_3_idx_plus:
                        fs = self.utr_3_idx_plus[itv[0]].find_gene(itv[1],itv[2])
            
                if strand == "-" or strand == ".":
                    if itv[0] in self.utr_3_idx_minus:
                        fs = fs + self.utr_3_idx_minus[itv[0]].find_gene(itv[1], itv[2])
                if len(fs) > 0:
                    genes = genes + fs
    
            except:
                raise
        genes = sorted(set(genes))
        return genes

    #
    def ovp_intron(self,itv_list,strand) :
        fs = []
        for (chrom,st,end) in itv_list :
            try :
                if strand == "+" or strand == ".":
                    if chrom in self.intron_idx_plus :
                        fs = self.intron_idx_plus[chrom].find_gene(st,end)
    
                if strand == "-" or strand == "." :
                    if chrom in self.intron_idx_minus :
                        fs = fs + self.intron_idx_minus[chrom].find_gene(st,end)
            except:
                raise
            if len(fs) > 0 :
                return True
            
        return False

    def ovp_intergenicUp(self,itv_list,strand) :
        fs = []
        for (chrom,st,end) in itv_list :
            try :
                if strand == "+" or strand == ".":
                    if chrom in self.intergenic_up_1kb_plus :
                        fs = self.intergenic_up_1kb_plus[chrom].find_gene(st,end)
            
                if strand == "-" or strand == "." :
                    if chrom in self.intergenic_up_1kb_minus :
                        fs = fs + self.intergenic_up_1kb_minus[chrom].find_gene(st,end)
            except:
                raise
            if len(fs) > 0 :
                return True
            
        return False

    def ovp_intergenicDown (self,itv_list,strand) :
        fs = []
        for (chrom,st,end) in itv_list :
            try :
                if strand == "+" or strand == ".":
                    if chrom in self.intergenic_up_1kb_plus :
                        fs = self.intergenic_up_1kb_plus[chrom].find_gene(st,end)
            
                if strand == "-" or strand == "." :
                    if chrom in self.intergenic_up_1kb_minus :
                        fs = fs + self.intergenic_up_1kb_minus[chrom].find_gene(st,end)
            except:
                raise
            if len(fs) > 0 :
                return True
            
        return False


    def getFeatures(self) :
        return self.features

    def Gene_annotation(self,itv_list,strand):
        genes = []
        fs = []
        cds = False
        for itv in itv_list :
            try:
                if strand == "+" or strand == ".":
                    if itv[0] in self.cds_exon_idx_plus :
                        fs = self.cds_exon_idx_plus[itv[0]].find_gene(itv[1],itv[2])

                if strand == "-" or strand == ".":
                    if itv[0] in self.cds_exon_idx_minus:
                        fs += self.cds_exon_idx_minus[itv[0]].find_gene(itv[1], itv[2])

                if len(fs) > 0:
                        genes = genes + fs
        
            except:
                    raise

        genes = sorted(set(genes))
        
        if len(genes) > 0:
            cds = True
            return (genes,cds,False,False,False,False,False)
        #sys.stderr.write("not cds exon\n")
        genes = self.ovp_utr5(itv_list,strand)
        if len(genes) > 0 :
            utr5 = True
            return (genes,False,utr5,False,False,False,False)
        #sys.stderr.write("not 5 utr\n")
        
        genes = self.ovp_utr3(itv_list,strand)
        if len(genes) > 0 :
            utr3 = True
            return (genes,False,False,utr3,False,False,False)
                
        #sys.stderr.write("not 3 utr\n")
        intron = self.ovp_intron(itv_list,strand)
        if intron :
            return (genes,False,False,False,intron,False,False)
        
        #sys.stderr.write("not intron\n")
        up1k = self.ovp_intergenicUp(itv_list,strand)
        if up1k :
            return (genes,False,False,False,False,up1k,False)
        #sys.stderr.write("not up1k\n")
        dn1k = self.ovp_intergenicDown(itv_list,strand)
        if dn1k :
            return (genes,False,False,False,False,False,dn1k)
        return(genes,False,False,False,False,False,False)







