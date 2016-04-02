CXX = g++
CXXFLAGS = -fopenmp -O3 -mtune=core2 -march=core2 -m64 -msse4.2 -g -ggdb -std=c++11
LDFLAGS = -fopenmp

OBJS = bench_cXc.o complexField.o basicmath.o

all: bench_cXc

cXc: $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o bench_cXc

complexField.o: ../complexField.cpp ../complexField.h
	$(CXX) $(CXXFLAGS) -c ../complexField.cpp

basicmath.o: ../basicmath.cpp ../basicmath.h
	$(CXX) $(CXXFLAGS) -c ../basicmath.cpp

cXc.o: bench_cXc.cpp ../basicmath.h
	$(CXX) $(CXXFLAGS) -c bench_cXc.cpp

clean:
	rm bench_cXc
	rm *.o

.PHONY: all clean
