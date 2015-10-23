# distutils: language = c++
# distutils: sources = GeneFeatures.cpp IntervalTree.cpp
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.pair cimport pair
#from libcpp.map cimport map

cdef extern from "IntervalTree.h" : #namespace "shapes":
    cdef cppclass IntervalTree
cdef extern from "GeneFeatures.h" :
#ctypedef struct chr_ITV :
#        string chrom
#        long start
#        long end
#
#ctypedef map[string, vector[string] ] gene_exon_Dict
#ctypedef map[string,map[int, IntervalTree ] ] chrom_itvTree_Dict

    cdef cppclass GeneFeatures:
        #vector[string] features
        
        #gene_exon_Dict plus_gene_exons
        #gene_exon_Dict minus_gene_exons
        
        #chrom_itvTree_Dict cds_exon_idx_plus
        #chrom_itvTree_Dict cds_exon_idx_minus
        
        GeneFeatures(string GTFfilename,string id_attribute)
        
        vector[pair[long,long] ] find_exons(string chr ,long st,long end,string strand,string gene)
        vector[string] getFeatures()
        vector[string] get_exon_pair(string g);
        #vector[string] Gene_annotation(vector[chr_ITV] itvlist, string strand)
        vector[string] Gene_annotation(string chrom, vector[pair[long,long]] itvlist, string strand)

cdef extern from "rRNA.h" :
    cdef cppclass rRNA :
        rRNA(string rRNAfilename)
        #bint is_rRNA(vector[chr_ITV] itvlist, string strand)
        bint is_rRNA(string chrom,vector[pair[long,long]] itvlist, string strand)

cdef class PyGeneFeatures:
    cdef GeneFeatures *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self,  GTFfilename,  id_attribute):
        self.thisptr = new GeneFeatures(GTFfilename, id_attribute)
    
    def __dealloc__(self):
        del self.thisptr
    
    def find_exons(self, chrom,st,end,strand,gene) :
        return self.thisptr.find_exons(chrom, st,end, strand, gene)
        
    def getFeatures(self) :
        return self.thisptr.getFeatures()
    def get_exon_pair(self, gene) :
        return self.thisptr.get_exon_pair(gene)

    def Gene_annotation(self, chrom,itv_list, string strand):
        return self.thisptr.Gene_annotation(chrom,itv_list,strand)



cdef class PyrRNA :
    cdef rRNA *thisptr
    def __cinit__(self,rRNAfilename):
        self.thisptr = new rRNA(rRNAfilename)

    def __dealloc__(self):
        del self.thisptr

    def is_rRNA(self, chrom,itv_list, strand) :
        return self.thisptr.is_rRNA(chrom, itv_list,strand)






