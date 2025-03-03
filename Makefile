CXX = g++
CXXFLAGS = -L/opt/homebrew/lib/ -I/opt/homebrew/include/ -std=c++20
LDFLAGS = -lCatch2 -lCatch2Main

all: test bench

test: test.cpp $(wildcard *.hpp)
	$(CXX) $(CXXFLAGS) test.cpp -o test $(LDFLAGS)

bench: bench.cpp discrete_laplacian_distribution.hpp canonne.hpp
	$(CXX) $(CXXFLAGS) bench.cpp -O3 -o bench $(LDFLAGS)

clean:
	rm test bench
