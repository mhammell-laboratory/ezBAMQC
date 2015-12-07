# Makefile for BAMqc, utilities for the Sequence Alignment/Map format.
#
#    Copyright (C) 2015 Bioinformatics Shared Resource, CSHL.
#    Portions copyright (C) 2015 Cold Spring Harbor Laboratory.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

CC       = g++
CPPFLAGS = $(DFLAGS) $(INCLUDES)
CFLAGS   = -g -fpermissive -Wall -O9 -O3 -std=c++11 -fPIC 
LDFLAGS  = -O9 -fpermissive
LDLIBS   =
DFLAGS=     -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_CURSES_LIB=1
LOBJS=      src/GeneFeatures.o src/rRNA.o src/IntervalTree.o src/InnerDist_prof.o \
            src/Results.o src/Mappability.o src/Coverage_prof.o src/parseBAM.o

INCLUDES=   -I./include -I$(HTSDIR)
LIBCURSES=  -lcurses # -lXCurses

prefix      = /usr/local
exec_prefix = $(prefix)
bindir      = $(exec_prefix)/bin
mandir      = $(prefix)/share/man
man1dir     = $(mandir)/man1

MKDIR_P = mkdir -p
#INSTALL = install -p
#INSTALL_PROGRAM = $(INSTALL)
#INSTALL_DATA    = $(INSTALL) -m 644
#INSTALL_DIR     = $(MKDIR_P) -m 755


PROGRAMS = libBAMqc.so


all: $(PROGRAMS)


# Adjust $(HTSDIR) to point to your top-level htslib directory
HTSDIR = ./htslib
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
BGZIP  = $(HTSDIR)/bgzip


PACKAGE_VERSION = 0.5


.SUFFIXES: .cpp .o

.cpp.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<


#lib:libbam.a

#libbam.a:$(LOBJS)
#	$(AR) -csru $@ $(LOBJS)

libBAMqc.so: $(LOBJS) $(HTSLIB)
	$(CC) -pthread $(LDFLAGS) -o $@ $(AOBJS) $(HTSLIB) $(LDLIBS) $(LIBCURSES) -lm -lz

libBAMqc.so: $(LOBJS) $(HTSLIB)
#	$(CC) -shared -Wl,-soname, libBAMqc.so.$(PACKAGE_VERSION) -lpthread $(LDFLAGS) -o $@ $(LOBJS) $(HTSLIB) $(LDLIBS) -lz -lm
	$(CC) -shared  -lpthread $(LDFLAGS) -o $@ $(LOBJS) $(HTSLIB) $(LDLIBS) -lz -lm
#	ln -sf $@ libBAMqc.so.$(PACKAGE_VERSION)

Constants_h = include/Constants.h
IntervalTree_h = include/IntervalTree.h $(Constants_h)
GeneFeatures_h = include/GeneFeatures.h $(Constants_h)
rRNA_h = include/rRNA.h $(IntervalTree_h) $(GeneFeatures_h)
Results_h = include/Resualts.h
Mappability_h = include/Mappability.h $(Constants_h)
InnerDist_prof_h = include/InnerDist_prof.h $(GeneFeatures_h)
Coverage_prof_h = include/Coverage_prof.h $(GeneFeatures_h)
parseBAM_h = include/parseBAM.h


IntervalTree.o: src/IntervalTree.cpp $(IntervalTree_h)
GeneFeatures.o: src/GeneFeatures.cpp $(GeneFeatures_h)
rRNA.o: src/rRNA.cpp $(rRNA_h)
Results.o: src/Results.cpp $(Results_h)
Mappability.o: src/Mappability.cpp $(Mappability_h) $(htslib_sam_h)
InnerDist_prof.o: src/InnerDist_prof.cpp $(InnerDist_prof_h) $(htslib_sam_h)
Coverage_prof.o: src/Coverage_prof.cpp $(Coverage_prof_h)
parseBAM.o: src/parseBAM.cpp $(parseBAM_h) $(htslib_sam_h) $(GeneFeatures_h) $(rRNA_h) $(Mappability_h) $(Coverage_prof_h) $(InnerDist_prof_h)


#install: $(PROGRAMS) $(BUILT_MISC_PROGRAMS)
#	$(INSTALL_DIR) $(DESTDIR)$(bindir) $(DESTDIR)$(man1dir)
#	$(INSTALL_PROGRAM) $(PROGRAMS) $(MISC_PROGRAMS) $(DESTDIR)$(bindir)
#	$(INSTALL_DATA) samtools.1 $(DESTDIR)$(man1dir)


mostlyclean:
	-rm -f src/*.o

clean: mostlyclean
	-rm -f $(PROGRAMS)

distclean: clean
	-rm -f TAGS

clean-all: clean


tags:
	ctags -f TAGS *.[ch] misc/*.[ch]


force:


.PHONY: all clean clean-all distclean force 
.PHONY: mostlyclean tags
