#
# PMTCalib makefile
# author : Leonidas N. Kalousis
# e-mail : leonidas.kalousis@gmail.com
#

dict = PMTCalibDict
lib = lib/libPMTCalib.so
rootmap = $(lib:.so=.rootmap)

srcs = $(wildcard src/*.cc)
head = $(wildcard src/*.h)
objs = $(srcs:.cc=.o)

cxx = g++
cxxflags = -g -W -O -Wall -Wno-deprecated -Werror -fPIC -std=c++1y

incflags = -I.
incflags += -I$(ROOTSYS)/include/ -I/$(GSL)/include/gsl/ -I/$(PMTCALIB)/src/ 

so = g++
soflags = -g -shared -fPIC
libs = $(ROOTLIBS) -lMinuit -lMinuit2 -lGeom -lXMLIO -lfftw3 -L/$(GSL)/lib -lgsl -lgslcblas 

all	: start $(dict).cc $(objs) $(lib) end

start	: 
	@echo ''		
	@echo ' * PMTCalib make ... '
	@echo ''
	@rm -f ./#* ./*~ ./*.*~	
	@rm -f ./src/#* ./src/*~ ./src/*.*~
	@rm -f ./mac/#* ./mac/*~ ./mac/*.*~
	@mkdir -p lib

$(dict).cc : 
	@rootcling -f $(dict).cc -s $(lib) -rml $(lib) -rmf $(rootmap) $(incflags) -c $(head) LinkDef.h 
	@echo ' * Building ( dic ) :' $@
	@echo ''

%.o	: %.cc	
	@$(cxx) $(cxxflags) $(incflags) -c -o $@ $<
	@echo ' * Building ( obj ) :' $@
	@echo ''

$(lib) 	: $(objs) $(dict).o
	@$(so) $(soflags) $(libs) -o $@ $(objs) $(dict).o
	@echo ' * Building ( lib ) :' $@
	@echo ''

end	:
	@echo ' * PMTCalib done !'
	@echo ''

clean	:	
	@echo ''	
	@echo ' * PMTCalib clean ...'
	@echo ''
	@rm -f ./#* ./*~ ./*.*~	
	@rm -f ./src/#* ./src/*~ ./src/*.*~
	@rm -f ./mac/#* ./mac/*~ ./mac/*.*~
	@rm -f $(dict)*.cc
	@rm -f $(dict)*.o
	@rm -f ./src/*.o
	@rm -f ./$(lib)
	@rm -f ./$(lib:.so=_rdict.pcm)
	@rm -f ./$(lib:.so=.rootmap)

fresh 	: clean all 
