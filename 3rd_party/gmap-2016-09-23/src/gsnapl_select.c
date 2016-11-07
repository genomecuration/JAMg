#include "bool.h"
#include "cpuid.h"
#include <errno.h>
#include <sys/stat.h>		/* For fstat */
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <libgen.h>
#include <unistd.h>

static bool
file_exists_p (char *filename) {
  struct stat sb;

  if (stat(filename,&sb) == 0) {
    return true;
  } else {
    return false;
  }
}

int
main (int argc, char *argv[]) {
  bool sse2_support_p;
  bool ssse3_support_p;
  bool sse4_1_support_p;
  bool sse4_2_support_p;
  bool avx2_support_p;
  char *dir, **new_argv;
  int i;
  int rc;

  CPUID_support(&sse2_support_p,&ssse3_support_p,&sse4_1_support_p,&sse4_2_support_p,&avx2_support_p);

  new_argv = (char **) malloc((argc + 1) * sizeof(char *));
  for (i = 1; i < argc; i++) {
    new_argv[i] = argv[i];
  }
  new_argv[argc] = (char *) NULL;

  if (index(argv[0],'/') == NULL) {
    /* Depend on path */
    /* Cannot use file_exists_p, since it won't search PATH */

    if (avx2_support_p == true) {
      new_argv[0] = "gsnapl.avx2";
      if ((rc = execvp(new_argv[0],new_argv)) == -1 && errno == ENOENT) {
	fprintf(stderr,"Note: %s does not exist.  For faster speed, may want to compile package on an AVX2 machine\n",new_argv[0]);
      } else {
	free(new_argv);
	return rc;
      }
    }

    if (sse4_2_support_p == true) {
      new_argv[0] = "gsnapl.sse42";
      if ((rc = execvp(new_argv[0],new_argv)) == -1 && errno == ENOENT) {
	fprintf(stderr,"Note: %s does not exist.  For faster speed, may want to compile package on an SSE4.2 machine\n",new_argv[0]);
      } else {
	free(new_argv);
	return rc;
      }
    }

    if (sse4_1_support_p == true) {
      new_argv[0] = "gsnapl.sse41";
      if ((rc = execvp(new_argv[0],new_argv)) == -1 && errno == ENOENT) {
	fprintf(stderr,"Note: %s does not exist.  For faster speed, may want to compile package on an SSE4.1 machine\n",new_argv[0]);
      } else {
	free(new_argv);
	return rc;
      }
    }

    if (ssse3_support_p == true) {
      new_argv[0] = "gsnapl.ssse3";
      if ((rc = execvp(new_argv[0],new_argv)) == -1 && errno == ENOENT) {
	fprintf(stderr,"Note: %s does not exist.  For faster speed, may want to compile package on an SSSE3 machine\n",new_argv[0]);
      } else {
	free(new_argv);
	return rc;
      }
    }

    if (sse2_support_p == true) {
      new_argv[0] = "gsnapl.sse2";
      if ((rc = execvp(new_argv[0],new_argv)) == -1 && errno == ENOENT) {
	fprintf(stderr,"Note: %s does not exist.  For faster speed, may want to compile package on an SSE2 machine\n",new_argv[0]);
      } else {
	free(new_argv);
	return rc;
      }
    }

    if (true) {
      new_argv[0] = "gsnapl.nosimd";
      rc = execvp(new_argv[0],new_argv);
      free(new_argv);
      return rc;
    }

  } else {
    dir = dirname(argv[0]);

    if (avx2_support_p == true) {
      new_argv[0] = (char *) malloc((strlen(dir) + strlen("/") + strlen("gsnapl.avx2") + 1) * sizeof(char));
      sprintf(new_argv[0],"%s/gsnapl.avx2",dir);
      if (file_exists_p(new_argv[0]) == false) {
	fprintf(stderr,"Note: %s does not exist.  For faster speed, may want to compile package on an AVX2 machine\n",new_argv[0]);
	free(new_argv[0]);
      } else {
	rc = execvp(new_argv[0],new_argv);
	free(new_argv[0]);
	free(new_argv);
	return rc;
      }
    }

    if (sse4_2_support_p == true) {
      new_argv[0] = (char *) malloc((strlen(dir) + strlen("/") + strlen("gsnapl.sse42") + 1) * sizeof(char));
      sprintf(new_argv[0],"%s/gsnapl.sse42",dir);
      if (file_exists_p(new_argv[0]) == false) {
	fprintf(stderr,"Note: %s does not exist.  For faster speed, may want to compile package on an SSE4.2 machine\n",new_argv[0]);
	free(new_argv[0]);
      } else {
	rc = execvp(new_argv[0],new_argv);
	free(new_argv[0]);
	free(new_argv);
	return rc;
      }
    }

    if (sse4_1_support_p == true) {
      new_argv[0] = (char *) malloc((strlen(dir) + strlen("/") + strlen("gsnapl.sse41") + 1) * sizeof(char));
      sprintf(new_argv[0],"%s/gsnapl.sse41",dir);
      if (file_exists_p(new_argv[0]) == false) {
	fprintf(stderr,"Note: %s does not exist.  For faster speed, may want to compile package on an SSE4.1 machine\n",new_argv[0]);
	free(new_argv[0]);
      } else {
	rc = execvp(new_argv[0],new_argv);
	free(new_argv[0]);
	free(new_argv);
	return rc;
      }
    }

    if (ssse3_support_p == true) {
      new_argv[0] = (char *) malloc((strlen(dir) + strlen("/") + strlen("gsnapl.ssse3") + 1) * sizeof(char));
      sprintf(new_argv[0],"%s/gsnapl.ssse3",dir);
      if (file_exists_p(new_argv[0]) == false) {
	fprintf(stderr,"Note: %s does not exist.  For faster speed, may want to compile package on an SSSE3 machine\n",new_argv[0]);
	free(new_argv[0]);
      } else {
	rc = execvp(new_argv[0],new_argv);
	free(new_argv[0]);
	free(new_argv);
	return rc;
      }
    }

    if (sse2_support_p == true) {
      new_argv[0] = (char *) malloc((strlen(dir) + strlen("/") + strlen("gsnapl.sse2") + 1) * sizeof(char));
      sprintf(new_argv[0],"%s/gsnapl.sse2",dir);
      if (file_exists_p(new_argv[0]) == false) {
	fprintf(stderr,"Note: %s does not exist.  For faster speed, may want to compile package on an SSE2 machine\n",new_argv[0]);
	free(new_argv[0]);
      } else {
	rc = execvp(new_argv[0],new_argv);
	free(new_argv[0]);
	free(new_argv);
	return rc;
      }
    }

    if (true) {
      new_argv[0] = (char *) malloc((strlen(dir) + strlen("/") + strlen("gsnapl.nosimd") + 1) * sizeof(char));
      sprintf(new_argv[0],"%s/gsnapl.nosimd",dir);
      if (file_exists_p(new_argv[0]) == false) {
	fprintf(stderr,"Note: %s does not exist.  For faster speed, may want to compile package on an non-SIMD machine\n",new_argv[0]);
	free(new_argv[0]);
      } else {
	rc = execvp(new_argv[0],new_argv);
	free(new_argv[0]);
	free(new_argv);
	return rc;
      }
    }
  }

  fprintf(stderr,"Error: appropriate GSNAPL version not found\n");
  free(new_argv);
  return -1;
}

