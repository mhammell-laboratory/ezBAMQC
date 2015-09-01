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


class BED_Reader( ):

   """Parse a BED file

   Pass the constructor either a file name or an iterator of lines of a
   BED file. If a file name is specified, it may refer to a gzip compressed
   file.

   Yields tuple of (gene_id,chrom,strand,start position,end position,type)

   """

   def __init__( self, filename):
      self.line_no = None
      self.filename = filename

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
          items = line.split('\t')
          seqname = items[0]
          start = items[1]
          end = items[2]
          feature = items[3]
          strand = items[5]
          #( seqname, start, end, feature,score, strand,leftovers ) = line.split("\t")
          

          yield (seqname, strand, int(start), int(end), feature)

      lines.close()
      self.line_no = None

   def get_line_number_string( self ):
      if self.line_no is None:
            return "file %s closed" % self.filename

      else:
         return "line %d of file %s" % ( self.line_no, self.filename )


class rRNA:
    """index of rRNA annotations.
        """
    def __init__ (self):

        self.featureIdxs_plus = {}
        self.featureIdxs_minus = {}
        #self.features = []
    def set(self,rRNAfilename) :
        self.read_features(rRNAfilename)

    # Reading & processing annotation files
    def read_features(self,rRNA_filename) :

        #dict of dicts since the builtin type doesn't support it for some reason
        temp_plus = collections.defaultdict(dict)
        temp_minus = collections.defaultdict(dict)
        
        #sys.stderr.write("start to read rrNA bed file\n")
        # read count of features in BED file
        bed = BED_Reader(rRNA_filename)  # (id, seqname, strand, int(start),
        
        i = 0
        counts = 0
        
        #sys.stderr.write("start the loop\n")
        try:
            for (seqname, strand, start, end, feature) in bed:
                    if seqname is None :
                        continue
                
                    counts += 1
                    #sys.stderr.write("count = " + str(counts)+"\n")
                    
                    try:
                        if strand == "+"  :
                            temp_plus[seqname].append(Interval(feature,start,end))
                    except:
                        temp_plus[seqname] = [Interval(feature,start,end)]
    
                    try:
                        if strand == "-" :
                            temp_minus[seqname].append(Interval(feature,start,end))
                    except :
                        temp_minus[seqname] = [Interval(feature,start,end)]

                    #save gene id
                                  #if f[4] not in self.features :
                                  #self.features.append(f[4])

                    i += 1
                    if i % 100000 == 0 :
                        sys.stderr.write("%d BED lines processed.\n" % i)
        except:
            sys.stderr.write("Error occured in %s.\n" % bed.get_line_number_string())
            raise

        if counts == 0 :
            sys.stderr.write("Warning: No features of type found in rRNA bed file.\n" )

        #build interval trees
        for each_chrom in temp_plus.keys():
            self.featureIdxs_plus[each_chrom] = IntervalTree(temp_plus[each_chrom])


        for each_chrom in temp_minus.keys():
            self.featureIdxs_minus[each_chrom] = IntervalTree(temp_minus[each_chrom])



    def is_rRNA(self,itv_list,strand):
        rRNAs = []
        fs = []
        for itv in itv_list :
            try:
                if strand == "+" :
                    if itv[0] in self.featureIdxs_plus :
                        fs = self.featureIdxs_plus[itv[0]].find_gene(itv[1],itv[2])


                if strand == "-" :
                        if itv[0] in self.featureIdxs_minus:
                            fs = self.featureIdxs_minus[itv[0]].find_gene(itv[1], itv[2])


                if strand == "." :
                        if itv[0] in self.featureIdxs_minus:
                            fs = self.featureIdxs_minus[itv[0]].find_gene(itv[1], itv[2])

                        if itv[0] in self.featureIdxs_plus :
                            fs += self.featureIdxs_plus[itv[0]].find_gene(itv[1],itv[2])
 

                if len(fs) > 0:
                        return True

            except:
                    raise


        return False
