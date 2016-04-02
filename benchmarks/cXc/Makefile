top_srcdir = ../..

CXX = g++
CXXFLAGS = -fopenmp -O3 -mtune=core2 -march=core2 -m64 -msse4.2 -g -ggdb -std=c++11
LDFLAGS = -fopenmp

OBJS = bench_cXc.o complexField.o basicmath.o

all: bench_cXc

bench_cXc: $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o bench_cXc

complexField.o: $(top_srcdir)/complexField.cpp $(top_srcdir)/complexField.h
	$(CXX) $(CXXFLAGS) -c $(top_srcdir)/complexField.cpp

basicmath.o: $(top_srcdir)/basicmath.cpp $(top_srcdir)/basicmath.h
	$(CXX) $(CXXFLAGS) -c $(top_srcdir)/basicmath.cpp

bench_cXc.o: bench_cXc.cpp $(top_srcdir)/basicmath.h
	$(CXX) $(CXXFLAGS) -c bench_cXc.cpp

clean:
	rm bench_cXc
	rm *.o

.PHONY: all clean
