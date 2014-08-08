#include <stdio.h>
#include "hsw.c"
int main() {
    if(can_use_intel_core_4th_gen_features())
	printf("Yes\n");
    else
	printf("No\n");
    return 0;
}
