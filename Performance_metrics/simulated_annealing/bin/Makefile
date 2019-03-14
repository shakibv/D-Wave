#
#	Makefile
#
###
# 
#	Copyright (C) 2013 Matthias Troyer <troyer@phys.ethz.ch>
#			   Ilia Zintchenko <zintchenko@itp.phys.ethz.ch>
#

.PHONY: all single threaded clean

.DEFAULT: all

CXX = g++

CXXFLAGS = -Wall -ansi -pedantic -std=c++11 -O3 -funroll-loops -pipe

TARGETS = an_ms_r1_fi an_ms_r1_nf an_ms_r1_nf_v0 an_ms_r3_nf an_ss_ge_fi an_ss_ge_fi_vdeg an_ss_ge_nf_bp an_ss_ge_nf_bp_vdeg an_ss_ge_fi_bp_vdeg an_ss_rn_fi an_ss_rn_fi_vdeg

TARGETS_OMP = $(addsuffix _omp,$(TARGETS))

all: single threaded

single: $(TARGETS)

threaded: $(TARGETS_OMP)

clean:
	rm -f $(TARGETS) $(TARGETS_OMP)

$(TARGETS) : %: main2.cc %.h sched.h usage.h utils.h output.h bits.h lattice.h
	$(CXX) $(CXXFLAGS) -DALGORITHM=\"$@.h\" -o $@ $<

$(TARGETS_OMP) : %: main_omp2.cc $(%.h:_omp=) sched.h usage.h utils.h output.h bits.h lattice.h
	$(CXX) $(CXXFLAGS) -fopenmp -DALGORITHM=\"$(@:_omp=).h\" -o $@ $<
