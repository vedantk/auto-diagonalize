#include <cstdio>
#include <cstdlib>
#include <cstdint>

double dragon(int n) {
    double a = 1.0;
    double b = 2.0;
    double c = 4.0;
	for (int i=3; i <= n; ++i) {
        double tmp = a;
        a = (.00001 * b) + (.00002 * c);
        b = (tmp / 2.0) + b;
        c = (tmp / 3.0) + c;
	}
	return c;
}

int main(int argc, char** argv) {
	if (argc != 2) {
		puts("Usage: ./dragon <number>");
		return 1;
	}

	int n = atoi(argv[1]);
	printf("dragon(%d) = %f\n", n, dragon(n));
	return 0;
}
