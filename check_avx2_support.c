#include <stdio.h>
int main() {
    if(__builtin_cpu_supports("avx2"))
	printf("Yes\n");
    else
	printf("No\n");
    return 0;
}
