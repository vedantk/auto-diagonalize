CXX = clang++
CXXFLAGS = -Wall -Wextra -std=c++11 -g `llvm-config --cxxflags` \
	-I/usr/include/eigen3

diagonalize.so: diagonalize.o
	$(CXX) -shared -Wl,-soname,$@ -o $@ $^

diagonalize.o: diagonalize.cc
	$(CXX) -c -fPIC $(CXXFLAGS) $^ -o $@

clean:
	rm -f test/*.bc
	rm -f test/*.bc.diag
	rm -f test/*.out
	rm -f diagonalize.{o,so}
