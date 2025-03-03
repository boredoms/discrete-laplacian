CXX = g++
CXXFLAGS = -std=c++20
INCLUDES = -I/opt/homebrew/include/
LDFLAGS = -L/opt/homebrew/lib/ 
LDLIBS = -lCatch2 -lCatch2Main

all: test bench

test: test.cpp $(wildcard *.hpp)
	$(CXX) $(CXXFLAGS) $(INCLUDES) test.cpp -o test $(LDFLAGS) $(LDLIBS)

bench: bench.cpp discrete_laplacian_distribution.hpp canonne.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) bench.cpp -O3 -o bench $(LDFLAGS) $(LDLIBS)

clean:
	rm test bench
