CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wpedantic

# the homebrew includes are for Mac OS users, you can remove them if you're compiling the tests on Linux
INCLUDES = -I/opt/homebrew/include/ -Iinclude

LDFLAGS = -L/opt/homebrew/lib/ 
LDLIBS = -lCatch2 -lCatch2Main

.PHONY: all test bench clean

all: test bench

test: test/test.cpp test/test.hpp include/discrete_distributions/$(wildcard *.hpp)
	$(CXX) $(CXXFLAGS) $(INCLUDES) test/test.cpp -o build/test $(LDFLAGS) $(LDLIBS)

bench: test/bench.cpp include/discrete_distributions/$(wildcard *.hpp)
	$(CXX) $(CXXFLAGS) $(INCLUDES) test/bench.cpp -O3 -o build/bench $(LDFLAGS) $(LDLIBS)

clean:
	rm build/*
