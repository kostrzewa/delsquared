top_srcdir = ../..

CXX = g++
CXXFLAGS = -fopenmp -O3 -mtune=core2 -march=core2 -m64 -msse4.2 -g -ggdb -std=c++11
LDFLAGS = -fopenmp

OBJS = bench_su3mXsu3v.o su3mField.o su3vField.o basicmath.o

all: bench_su3mXsu3v

bench_su3mXsu3v: $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o bench_su3mXsu3v

su3mField.o: $(top_srcdir)/su3mField.cpp $(top_srcdir)/su3mField.h
	$(CXX) $(CXXFLAGS) -c $(top_srcdir)/su3mField.cpp

su3vField.o: $(top_srcdir)/su3vField.cpp $(top_srcdir)/su3vField.h
	$(CXX) $(CXXFLAGS) -c $(top_srcdir)/su3vField.cpp

basicmath.o: $(top_srcdir)/basicmath.cpp $(top_srcdir)/basicmath.h
	$(CXX) $(CXXFLAGS) -c $(top_srcdir)/basicmath.cpp

bench_su3Xsu3v.o: bench_su3mXsu3v.cpp $(top_srcdir)/basicmath.h $(top_srcdir)/su3vField.h $(top_srcdir)/su3mField.h
	$(CXX) $(CXXFLAGS) -c bench_su3mXsu3v.cpp

clean:
	rm bench_su3mXsu3v
	rm *.o

.PHONY: all clean
