
#
# makefile
#
# author: Leonidas N. Kalousis
#

# files

srcs = ana1.cc

objs = $(srcs:.cc=.o)

exe = $(srcs:.cc=)-exec

# libs

cxx = g++

cxxflags = -g -W -O -Wall -Wno-deprecated -fPIC -std=c++14

cxxflags += -I$(ROOTSYS)/include/ -I$(PMTCALIB)/src/ -I/usr/local/include/gsl/

rootlibs = $(shell root-config --glibs) -lMinuit -lMinuit2 -lfftw3 -lm -lUnfold -L/$(GSL)/ -lgsl -lgslcblas -L$(PMTCALIB)/lib -lPMTCalib

# targets

all	: start $(exe) end

start	: 
	@echo ''		
	@echo ' * make ...'
	@echo ''
	@rm -f ./#* ./*~ ./*.*~	

$(exe)	: $(objs) 
	@$(cxx) $(rootlibs) $(cxxflags) -o $@ $^
	@echo ' * building ( exe )  :' $@
	@echo ''

%.o	: %.cc	
	@$(cxx) $(cxxflags) -c -o $@ $<
	@echo ' * building ( objs ) :' $@
	@echo ''

end	:
	@echo ' * make done !'
	@echo ''

clean	:	
	@echo ''	
	@echo ' * make clean ...'
	@echo ''		
	@rm -f ./#* ./*~ ./*.*~ *.o *-exec

fresh	: clean all
