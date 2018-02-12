mydefs = -DNPT
incl = -I fp
objs = fpdim3.o fprandom.o fpvtfhandler.o mcmover.o spherocyl.o fpgrid.o simulation.o cosinator.o icmaker.o vtfsaver.o kdtree.o fpspace.o kdtreespherocyl.o vtfloader.o analyzer.o
links = -lgsl -lgslcblas -lm
CXXFLAGS = 
CCFLAGS = 


default: $(objs) testrun.o kdtest.o spcrun.o analyze.o
	$(CXX) $(CXXFLAGS) -o testrun $(objs) testrun.o $(links) $(incl)
	$(CXX) $(CXXFLAGS) -o spcrun $(objs) spcrun.o $(links) $(incl)
	$(CXX) $(CXXFLAGS) -o kdtest $(objs) kdtest.o $(links) $(incl)
	$(CXX) $(CXXFLAGS) -o analyze $(objs) analyze.o $(links) $(incl)
	
fpdim3.o: fp/fpdim3.cpp
	$(CXX) $(CXXFLAGS) $(mydefs) $(incl) -c fp/fpdim3.cpp

fpvtfhandler.o: fp/fpvtfhandler.cpp
	$(CXX) $(CXXFLAGS) $(mydefs) $(incl) -c fp/fpvtfhandler.cpp

fprandom.o: fp/fprandom.cpp
	$(CXX) $(CXXFLAGS) $(mydefs) $(incl) -c fp/fprandom.cpp

fpStatistics.o: fp/fpStatistics.cpp
	$(CXX) $(CXXFLAGS) $(mydefs) $(incl) -c fp/fpStatistics.cpp

testrun.o: testrun.cpp
	$(CXX) $(CXXFLAGS) $(mydefs) $(incl) -c testrun.cpp

spcrun.o: spcrun.cpp
	$(CXX) $(CXXFLAGS) $(mydefs) $(incl) -c spcrun.cpp

kdtest.o: kdtest.cpp
	$(CXX) $(CXXFLAGS) $(mydefs) $(incl) -c kdtest.cpp

mcmover.o: mcmover.cpp
	$(CXX) $(CXXFLAGS) $(mydefs) $(incl) -c mcmover.cpp

spherocyl.o: spherocyl.cpp
	$(CXX) $(CXXFLAGS) $(mydefs) $(incl) -c spherocyl.cpp

fpgrid.o: fpgrid.cpp
	$(CXX) $(CXXFLAGS) $(mydefs) $(incl) -c fpgrid.cpp

simulation.o: simulation.cpp
	$(CXX) $(CXXFLAGS) $(mydefs) $(incl) -c simulation.cpp

cosinator.o: cosinator.cpp
	$(CXX) $(CXXFLAGS) $(mydefs) $(incl) -c cosinator.cpp

icmaker.o: icmaker.cpp
	$(CXX) $(CXXFLAGS) $(mydefs) $(incl) -c icmaker.cpp

vtfsaver.o: vtfsaver.cpp
	$(CXX) $(CXXFLAGS) $(mydefs) $(incl) -c vtfsaver.cpp

vtfloader.o: vtfloader.cpp
	$(CXX) $(CXXFLAGS) $(mydefs) $(incl) -c vtfloader.cpp

kdtree.o: kdtree.cpp
	$(CXX) $(CXXFLAGS) $(mydefs) $(incl) -c kdtree.cpp

kdtreespherocyl.o: kdtreespherocyl.cpp
	$(CXX) $(CXXFLAGS) $(mydefs) $(incl) -c kdtreespherocyl.cpp

fpspace.o: fpspace.cpp
	$(CXX) $(CXXFLAGS) $(mydefs) $(incl) -c fpspace.cpp

analyzer.o: analyzer.cpp
	$(CXX) $(CXXFLAGS) $(mydefs) $(incl) -c analyzer.cpp

analyze.o: analyze.cpp
	$(CXX) $(CXXFLAGS) $(mydefs) $(incl) -c analyze.cpp

clean:
	rm -f analyze main testrun spcrun *.o fp/*.o

debug: CXXFLAGS += -Ddebug -g3 -g
debug: CCFLAGS += -Ddebug -g3 -g
debug: default

fast: CXXFLAGS += -O3
fast: CCFLAGS += -O3
fast: default

warning: CXXFLAGS += -Wall
warning: CCFLAGS += -Wall
warning: default

NPT: CXXFLAGS += -DNPT -g3 -g
NPT: CCFLAGS += -DNPT -g3 -g
NPT: default
