test: test.cpp dld.hpp
	g++ test.cpp -L/opt/homebrew/lib/ -I/opt/homebrew/include/ -std=c++20 -lCatch2 -lCatch2Main -o test

clean:
	rm test
