#include <cstdio>
#include <cstdlib>
#include <cstdint>

void foo(int n) {
	if (n % 100 == 0) {
		puts("Beep.");
	}
}

int dragon(int n) {
	int a = 8;
	int b = 4;
	int c = 2;
	for (int i=3; i <= n; ++i) {
		int ap = c;
		int bp = b + (a / 2);
		int cp = a;
		a = ap;
		b = bp;
		c = cp;
		foo(n);
	}
	return a;
}

int main(int argc, char** argv) {
	if (argc != 2) {
		puts("Usage: ./dragon <number>");
		return 1;
	}

	int n = atoi(argv[1]);
	printf("dragon(%d) = %d\n", n, dragon(n));
	return 0;
}
