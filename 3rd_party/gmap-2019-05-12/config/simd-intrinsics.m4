
# _mm_extract_epi64 is supposed to be covered by HAVE_SSE4_1, but fails on i386 machines


AC_DEFUN([ACX_SIMD_INTRINSICS], [
AC_LANG_SAVE
AC_LANG(C)

AC_MSG_CHECKING(for _mm_extract_epi64)
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([[#include <smmintrin.h>]],
                   [[__m128i a;]],
                   [[_mm_extract_epi64(a,0);]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_MM_EXTRACT_EPI64],[1],[Define to 1 if _mm_extract_epi64 intrinsic is available.])],
  [AC_MSG_RESULT(no)])

AC_MSG_CHECKING(for _mm_popcnt_u64)
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([[#include <smmintrin.h>]],
                   [[int a;]],
                   [[_mm_popcnt_u64(a);]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_MM_POPCNT_U64],[1],[Define to 1 if _mm_popcnt_u64 intrinsic is available.])],
  [AC_MSG_RESULT(no)])

AC_LANG_RESTORE
])


