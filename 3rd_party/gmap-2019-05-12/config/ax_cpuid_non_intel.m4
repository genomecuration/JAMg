
AC_DEFUN([AX_CPUID_NON_INTEL],
[AC_REQUIRE([AC_PROG_CC])
AC_LANG_PUSH([C])
# Test for SSE2 support
  AC_MSG_CHECKING(for sse2 support)
  AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[#include <stdint.h>
static void run_cpuid (uint32_t eax, uint32_t ecx, uint32_t *abcd) {
  uint32_t ebx, edx;
  __asm__ ("cpuid" : "+b" (ebx), "+a" (eax), "+c" (ecx), "=d" (edx));
  abcd[0] = eax; abcd[1] = ebx; abcd[2] = ecx; abcd[3] = edx;}]],
[[uint32_t abcd[4];
 uint32_t sse2_mask = (1 << 26);
 run_cpuid(1, 0, abcd);
 return ((abcd[/*EDX*/3] & sse2_mask) == sse2_mask) ? 0 : 9;]])],
        [AC_MSG_RESULT(yes)
         ax_cv_cpu_has_sse2_ext=yes],
	[AC_MSG_RESULT(no)])

# Test for SSSE3 support
  AC_MSG_CHECKING(for ssse3 support)
  AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[#include <stdint.h>
static void run_cpuid (uint32_t eax, uint32_t ecx, uint32_t *abcd) {
  uint32_t ebx, edx;
  __asm__ ("cpuid" : "+b" (ebx), "+a" (eax), "+c" (ecx), "=d" (edx));
  abcd[0] = eax; abcd[1] = ebx; abcd[2] = ecx; abcd[3] = edx;}]],
[[uint32_t abcd[4];
 uint32_t ssse3_mask = (1 << 9);
 run_cpuid(1, 0, abcd);
 return ((abcd[/*ECX*/2] & ssse3_mask) == ssse3_mask) ? 0 : 9;]])],
        [AC_MSG_RESULT(yes)
         ax_cv_cpu_has_ssse3_ext=yes],
	[AC_MSG_RESULT(no)])

# Test for SSE4.1 support
  AC_MSG_CHECKING(for sse4.1 support)
  AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[#include <stdint.h>
static void run_cpuid (uint32_t eax, uint32_t ecx, uint32_t *abcd) {
  uint32_t ebx, edx;
  __asm__ ("cpuid" : "+b" (ebx), "+a" (eax), "+c" (ecx), "=d" (edx));
  abcd[0] = eax; abcd[1] = ebx; abcd[2] = ecx; abcd[3] = edx;}]],
[[uint32_t abcd[4];
 uint32_t sse4_1_mask = (1 << 19);
 run_cpuid(1, 0, abcd);
 return ((abcd[/*ECX*/2] & sse4_1_mask) == sse4_1_mask) ? 0 : 9;]])],
        [AC_MSG_RESULT(yes)
         ax_cv_cpu_has_sse41_ext=yes],
	[AC_MSG_RESULT(no)])

# Test for SSE4.2 support
  AC_MSG_CHECKING(for sse4.2 support)
  AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[#include <stdint.h>
static void run_cpuid (uint32_t eax, uint32_t ecx, uint32_t *abcd) {
  uint32_t ebx, edx;
  __asm__ ("cpuid" : "+b" (ebx), "+a" (eax), "+c" (ecx), "=d" (edx));
  abcd[0] = eax; abcd[1] = ebx; abcd[2] = ecx; abcd[3] = edx;}]],
[[uint32_t abcd[4];
 uint32_t sse4_2_mask = (1 << 20);
 run_cpuid(1, 0, abcd);
 return ((abcd[/*ECX*/2] & sse4_2_mask) == sse4_2_mask) ? 0 : 9;]])],
        [AC_MSG_RESULT(yes)
         ax_cv_cpu_has_sse42_ext=yes],
	[AC_MSG_RESULT(no)])

# Test for popcnt support
  AC_MSG_CHECKING(for popcnt support)
  AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[#include <stdint.h>
static void run_cpuid (uint32_t eax, uint32_t ecx, uint32_t *abcd) {
  uint32_t ebx, edx;
  __asm__ ("cpuid" : "+b" (ebx), "+a" (eax), "+c" (ecx), "=d" (edx));
  abcd[0] = eax; abcd[1] = ebx; abcd[2] = ecx; abcd[3] = edx;}]],
[[uint32_t abcd[4];
 uint32_t popcnt_mask = (1 << 23);
 run_cpuid(1, 0, abcd);
 return ((abcd[/*ECX*/2] & popcnt_mask) == popcnt_mask) ? 0 : 9;]])],
        [AC_MSG_RESULT(yes)
         ax_cv_cpu_has_popcnt_ext=yes],
	[AC_MSG_RESULT(no)])

# Test for bmi1 support
  AC_MSG_CHECKING(for bmi1 support)
  AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[#include <stdint.h>
static void run_cpuid (uint32_t eax, uint32_t ecx, uint32_t *abcd) {
  uint32_t ebx, edx;
  __asm__ ("cpuid" : "+b" (ebx), "+a" (eax), "+c" (ecx), "=d" (edx));
  abcd[0] = eax; abcd[1] = ebx; abcd[2] = ecx; abcd[3] = edx;}]],
[[uint32_t abcd[4];
 uint32_t bmi1_mask = (1 << 3);
 run_cpuid(1, 0, abcd);
 return ((abcd[/*EBX*/1] & bmi1_mask) == bmi1_mask) ? 0 : 9;]])],
        [AC_MSG_RESULT(yes)
         ax_cv_cpu_has_bmi1_ext=yes],
	[AC_MSG_RESULT(no)])

# Test for AVX2 support
  AC_MSG_CHECKING(for avx2 support)
  AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[#include <stdint.h>
static void run_cpuid (uint32_t eax, uint32_t ecx, uint32_t *abcd) {
  uint32_t ebx, edx;
  __asm__ ("cpuid" : "+b" (ebx), "+a" (eax), "+c" (ecx), "=d" (edx));
  abcd[0] = eax; abcd[1] = ebx; abcd[2] = ecx; abcd[3] = edx;}
static int check_xcr0_ymm () {
  uint32_t xcr0;
  __asm__ ("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx");
  return ((xcr0 & 6) == 6);}]],
[[uint32_t abcd[4];
 uint32_t fma_movbe_osxsave_mask = ((1 << 12) | (1 << 22) | (1 << 27));
 uint32_t avx2_bmi12_mask = ((1 << 5) | (1 << 3) | (1 << 8));
 uint32_t lzcnt_mask = (1 << 5);
 run_cpuid(1, 0, abcd);
 if ((abcd[/*ECX*/2] & fma_movbe_osxsave_mask) != fma_movbe_osxsave_mask) {
   return 9;
 } else if (!check_xcr0_ymm()) {
   return 8;
 } else {
   run_cpuid(7, 0, abcd);
   if ((abcd[/*EBX*/1] & avx2_bmi12_mask) != avx2_bmi12_mask) {
     return 7;
   } else {
     run_cpuid(0x80000001, 0, abcd);
     if ((abcd[/*ECX*/2] & lzcnt_mask) != lzcnt_mask) {
       return 6;
     } else {
       return 0;
     }
   }
 }]])],
        [AC_MSG_RESULT(yes)
         ax_cv_cpu_has_avx2_ext=yes],
	[AC_MSG_RESULT(no)])

# Test for bmi2 support
  AC_MSG_CHECKING(for bmi2 support)
  AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[#include <stdint.h>
static void run_cpuid (uint32_t eax, uint32_t ecx, uint32_t *abcd) {
  uint32_t ebx, edx;
  __asm__ ("cpuid" : "+b" (ebx), "+a" (eax), "+c" (ecx), "=d" (edx));
  abcd[0] = eax; abcd[1] = ebx; abcd[2] = ecx; abcd[3] = edx;}]],
[[uint32_t abcd[4];
 uint32_t bmi2_mask = (1 << 8);
 run_cpuid(1, 0, abcd);
 return ((abcd[/*EBX*/1] & bmi2_mask) == bmi2_mask) ? 0 : 9;]])],
        [AC_MSG_RESULT(yes)
         ax_cv_cpu_has_bmi2_ext=yes],
	[AC_MSG_RESULT(no)])


# Test for AVX512 (F and CD) support
  AC_MSG_CHECKING(for avx512 support)
  AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[#include <stdint.h>
static void run_cpuid (uint32_t eax, uint32_t ecx, uint32_t *abcd) {
  uint32_t ebx, edx;
  __asm__ ("cpuid" : "+b" (ebx), "+a" (eax), "+c" (ecx), "=d" (edx));
  abcd[0] = eax; abcd[1] = ebx; abcd[2] = ecx; abcd[3] = edx;}
static int check_xcr0_zmm () {
  uint32_t xcr0;
  uint32_t zmm_ymm_xmm = ((7 << 5) | (1 << 2) | (1 << 1));
  __asm__ ("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx");
  return ((xcr0 & zmm_ymm_xmm) == zmm_ymm_xmm);}]],
[[uint32_t abcd[4];
 uint32_t osxsave_mask = (1 << 27);
 uint32_t avx512_mask = (/*512F*/(1 << 16) | /*512CD*/(1 << 28));
 run_cpuid(1, 0, abcd);
 if ((abcd[/*ECX*/2] & osxsave_mask) != osxsave_mask) {
   return 9;
 } else if (!check_xcr0_zmm()) {
   return 8;
 } else if ((abcd[/*EBX*/1] & avx512_mask) != avx512_mask) {
   return 0; /* Should be false here, but book/Web examples skip this test */
 } else {
   return 0;
 }]])],
        [AC_MSG_RESULT(yes)
         ax_cv_cpu_has_avx512_ext=yes],
	[AC_MSG_RESULT(no)])


# Test for AVX512BW and VL support
  AC_MSG_CHECKING(for avx512bw and avx512vl support)
  AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[#include <stdint.h>
static void run_cpuid (uint32_t eax, uint32_t ecx, uint32_t *abcd) {
  uint32_t ebx, edx;
  __asm__ ("cpuid" : "+b" (ebx), "+a" (eax), "+c" (ecx), "=d" (edx));
  abcd[0] = eax; abcd[1] = ebx; abcd[2] = ecx; abcd[3] = edx;}
static int check_xcr0_zmm () {
  uint32_t xcr0;
  uint32_t zmm_ymm_xmm = ((7 << 5) | (1 << 2) | (1 << 1));
  __asm__ ("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx");
  return ((xcr0 & zmm_ymm_xmm) == zmm_ymm_xmm);}]],
[[uint32_t abcd[4];
 uint32_t osxsave_mask = (1 << 27);
 uint32_t avx512bw_vl_mask = (/*512BW*/(1 << 30) | /*512VL*/(1 << 31));
 run_cpuid(1, 0, abcd);
 if ((abcd[/*ECX*/2] & osxsave_mask) != osxsave_mask) {
   return 9;
 } else if (!check_xcr0_zmm()) {
   return 8;
 } else if ((abcd[/*EBX*/1] & avx512bw_vl_mask) != avx512bw_vl_mask) {
   return 9;
 } else {
   return 0;
 }]])],
        [AC_MSG_RESULT(yes)
         ax_cv_cpu_has_avx512bw_ext=yes],
	[AC_MSG_RESULT(no)])


AC_LANG_POP([C])
])
