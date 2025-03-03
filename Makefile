test: test.cpp discrete_laplacian_distribution.hpp discrete_gaussian_distribution.hpp test.hpp canonne.hpp
	g++ test.cpp -L/opt/homebrew/lib/ -I/opt/homebrew/include/ -std=c++20 -lCatch2 -lCatch2Main -o test

bench: bench.cpp discrete_laplacian_distribution.hpp canonne.hpp
	g++ bench.cpp -L/opt/homebrew/lib/ -I/opt/homebrew/include/ -std=c++20 -O3 -lCatch2 -lCatch2Main -o bench

clean:
	rm test bench
