#include <cstdio>

void fail_nested() {
	for (int i=0; i < 10; ++i) {
		for (int j=0; j < 10; ++j) {
			printf("%d", i + j);
		}
	}
}

void fail_linear() {
	int a = 1;
	int b = 2;
	for (int i=0; i < 10; ++i) {
		a = a * a;
		b = b / b;
		printf("%d", a + b);
	}
}

void fail_nostate() {
	for (int i=0; i < 10; ++i) {
		puts(".");
	}
}

int main() {
	fail_nested();
	fail_linear();
	fail_nostate();
	return 0;
}
