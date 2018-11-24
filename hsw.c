#include <stdint.h>


void run_cpuid(uint32_t eax, uint32_t ecx, uint32_t* abcd) {
    uint32_t ebx, edx;

    __asm__("cpuid" : "+b" (ebx), "+a" (eax), "+c" (ecx), "=d" (edx));
    abcd[0] = eax;
    abcd[1] = ebx;
    abcd[2] = ecx;
    abcd[3] = edx;
}

int check_xcr0_ymm() {
    uint32_t xcr0;
    __asm__("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx");
    return ((xcr0 & 6) == 6); /* checking if xmm and ymm state are enabled in XCR0 */
}

int check_4th_gen_intel_core_features() {
    uint32_t abcd[4];
    uint32_t fma_movbe_osxsave_mask = ((1 << 12) | (1 << 22) | (1 << 27));
    uint32_t avx2_bmi12_mask = (1 << 5) | (1 << 3) | (1 << 8);
    /* CPUID.(EAX=01H, ECX=0H):ECX.FMA[bit 12]==1 &&
CPUID.(EAX=01H, ECX=0H):ECX.MOVBE[bit 22]==1 &&
CPUID.(EAX=01H, ECX=0H):ECX.OSXSAVE[bit 27]==1 */
    run_cpuid(1, 0, abcd);
    if ((abcd[2] & fma_movbe_osxsave_mask) != fma_movbe_osxsave_mask)
        return 0;
    if (!check_xcr0_ymm())
        return 0;
    /* CPUID.(EAX=07H, ECX=0H):EBX.AVX2[bit 5]==1 &&
CPUID.(EAX=07H, ECX=0H):EBX.BMI1[bit 3]==1 &&
CPUID.(EAX=07H, ECX=0H):EBX.BMI2[bit 8]==1 */
    run_cpuid(7, 0, abcd);
    if ((abcd[1] & avx2_bmi12_mask) != avx2_bmi12_mask)
        return 0;
    /* CPUID.(EAX=80000001H):ECX.LZCNT[bit 5]==1 */
    run_cpuid(0x80000001, 0, abcd);
    if ((abcd[2] & (1 << 5)) == 0)
        return 0;
    return 1;
}


static int can_use_intel_core_4th_gen_features() {
    static int the_4th_gen_features_available = -1;
    /* test is performed once */
    if (the_4th_gen_features_available < 0)
        the_4th_gen_features_available = check_4th_gen_intel_core_features();
    return the_4th_gen_features_available;
}

int is_cpuid_ecx_bit_set(int eax, int bitidx) {
    int ecx = 0, edx = 0, ebx = 0;
    __asm__("cpuid"
            : "=b" (ebx),
            "=c" (ecx),
            "=d" (edx)
            : "a" (eax)
            );
    return (((ecx >> bitidx)&1) == 1);
}

int is_sse42_supported() {
    return is_cpuid_ecx_bit_set(1, 20);
}

/*   
#include <stdio.h>
int main(int argc, char** argv)
{
  if ( can_use_intel_core_4th_gen_features() )
    printf("This CPU supports ISA extensions introduced in Haswell\n");
  else
    printf("This CPU does not support all ISA extensions introduced in Haswell\n");
  return 1;
}
 */

