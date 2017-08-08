
AC_DEFUN([ACX_BUILTIN_POPCOUNT], [
AC_REQUIRE([AC_CANONICAL_HOST])
AC_LANG_SAVE
AC_LANG(C)

CFLAGS_ORIG=$CFLAGS

if test x"$ax_cv_c_compiler_vendor" = xintel; then
  TEST_CFLAGS="$CFLAGS_ORIG"
  POPCNT_CFLAGS=""
else
  TEST_CFLAGS="$CFLAGS_ORIG -mpopcnt"
  AX_CHECK_COMPILE_FLAG([$TEST_CFLAGS], [ax_cv_compile_popcnt_ext=yes], [ax_cv_ext_compile_problem=yes])
  if test x"$ax_cv_compile_popcnt_ext" != xyes; then
    AC_MSG_WARN([Your compiler does not support -mpopcnt])
    POPCNT_CFLAGS=""
  else
    CFLAGS="$CFLAGS_ORIG -mpopcnt"
    AC_LINK_IFELSE([AC_LANG_PROGRAM([])],
                   [ax_cv_link_popcnt=yes],
   	           [ax_cv_ext_linker_problem=yes])
    if test x"$ax_cv_link_popcnt" != xyes; then
      AC_MSG_WARN([Your compiler supports -mpopcnt but not your linker.  Can you try another linker or update yours?])
      POPCNT_CFLAGS=""
    else
      POPCNT_CFLAGS="-mpopcnt"
    fi
  fi
fi


# Test for __builtin functions with or without the -mpopcnt compiler flag
CFLAGS="$CFLAGS_ORIG $POPCNT_CFLAGS"

AC_MSG_CHECKING(for __builtin_popcount)
AC_RUN_IFELSE(
  [AC_LANG_PROGRAM([[]],
                   [[return (__builtin_popcount(0xffffffffu) == 32) ? 0 : 9;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_BUILTIN_POPCOUNT],[1],[Define to 1 if __builtin_popcount works.])],
  [AC_MSG_RESULT(no)])

AC_MSG_CHECKING(for __builtin_clz)
AC_RUN_IFELSE(
  [AC_LANG_PROGRAM([[]],
                   [[return (__builtin_clz(0x1u) == 31) ? 0 : 9;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_BUILTIN_CLZ],[1],[Define to 1 if __builtin_clz works.])],
  [AC_MSG_RESULT(no)])

AC_MSG_CHECKING(for __builtin_ctz)
AC_RUN_IFELSE(
  [AC_LANG_PROGRAM([[]],
                   [[return (__builtin_ctz(0x80000000u) == 31) ? 0 : 9;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_BUILTIN_CTZ],[1],[Define to 1 if __builtin_ctz works.])],
  [AC_MSG_RESULT(no)])

CFLAGS=$CFLAGS_ORIG

AC_SUBST(POPCNT_CFLAGS)

AC_LANG_RESTORE
])


