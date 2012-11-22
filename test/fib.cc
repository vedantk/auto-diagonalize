#include <cstdio>
#include <cstdlib>
#include <cstdint>

double fib(int n) {
	double a = 1;
	double b = 1;
	for (int i=2; i <= n; ++i) {
		double tmp = a;
		a = b;
		b = tmp + b;
	}
	return b;
}

int main(int argc, char** argv) {
    if (argc != 2) {
    	puts("Usage: ./fib <number>");
    	return 1;
    }
    
    int n = atoi(argv[1]);
    for (int i=0; i < 10000; ++i) {
        printf("fib(%d) = %lf\n", n, fib(n));
    }
    
    return 0;
}
