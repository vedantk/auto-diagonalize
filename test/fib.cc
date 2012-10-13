#include <cstdio>
#include <cstdlib>
#include <cstdint>

int fib(int n) {
	int a = 1;
	int b = 1;
	for (int i=2; i <= n; ++i) {
		int tmp = a;
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
	printf("fib(%d) = %d\n", n, fib(n));
	return 0;
}
