/* pegas.c    2019-10-03 */

/* Copyright 2015-2019 Emmanuel Paradis */

/* This file is part of the R-package `pegas'. */
/* See the file ../DESCRIPTION for licensing issues. */

//#include <R.h>
#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

/* returns 1 if both bases are different surely, 0 otherwise */
#define DifferentBase(a, b) (a & b) < 16

void haplotype_DNAbin(unsigned char *x, int *n, int *s, int *haplo)
{
  int i1, i2, s1, s2, flag, k;

  i1 = 0;
  while (i1 < *n - 1) {
    if (!haplo[i1]) {
      i2 = i1 + 1;
      flag = 1; /* initially the two sequences are considered identical */
      while (i2 < *n && flag > 0 ) { //modified to avoid break
        if (!haplo[i2]) {
          s1 = i1;
          s2 = i2;
          k = 0;
          //flag = 1; /* initially the two sequences are considered identical */
          while (k < *s) {
            //added below
            if (x[s1] != x[s2]){
              if (x[s1] == 0x04 || x[s2] == 0x04){
                flag = 0;
              }
              else if (x[s1] > 0x07 && x[s2] > 0x07){
                flag = 0;
              }
            }
            //end insertion, begin deletion
            //if (x[s1] != x[s2] && x[s1] > 0x07 && x[s2] > 0x07) { /* fix the fact that a and b could be gaps (-) or completely missing (?) */
            //    if (DifferentBase(x[s1], x[s2])) {
            //	flag = 0;
            //	break;
            //    }
            //}
            s1 += *n;
            s2 += *n;
            k++;
          }
          if (flag) haplo[i2] = i1 + 1;
        }
        i2++;
      }
    }
    i1++;
  }
}

void distDNA_pegas(unsigned char *x, int *n, int *s, double *d)
{
    int i1, i2, s1, s2, target, Nd;
    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
	for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1 += *n, s2 += *n) {
		if (x[s1] <= 0x07 || x[s2] <= 0x07) continue;
		if (DifferentBase(x[s1], x[s2])) Nd++;
	    }
	    d[target] = ((double) Nd);
	    target++;
	}
    }
}

SEXP unique_haplotype_loci(SEXP x, SEXP NROW, SEXP NCOL)
{
    SEXP res;
    int nr, nc, i, j, k, *H, flag;

    PROTECT(x = coerceVector(x, STRSXP));
    PROTECT(NROW = coerceVector(NROW, INTSXP));
    PROTECT(NCOL = coerceVector(NCOL, INTSXP));

    nr = INTEGER(NROW)[0];
    nc = INTEGER(NCOL)[0];

    PROTECT(res = allocVector(INTSXP, nc));
    H = INTEGER(res);
    memset(H, 0, nc * sizeof(int));

    j = 0;
    while (j < nc - 1) {
	if (!H[j]) {
	    k = j + 1;
	    while (k < nc) {
		if (!H[k]) {
		    i = 0;
		    flag = 1; /* initially the two haplotypes are considered identical */
		    while (i < nr) {
			if (strcmp(CHAR(STRING_ELT(x, i + j*nr)), CHAR(STRING_ELT(x, i + k*nr)))) { /* the strings are different */
			    flag = 0;
			    break;
			}
			i++;
		    }
		    if (flag) H[k] = j + 1;
		}
		k++;
	    }
	}
	j++;
    }

    UNPROTECT(4);
    return res;
}

static R_CMethodDef C_entries[] = {
    {"haplotype_DNAbin", (DL_FUNC) &haplotype_DNAbin, 4},
    {"distDNA_pegas", (DL_FUNC) &distDNA_pegas, 4},
    {NULL, NULL, 0}
};

SEXP read_bin_pegas(SEXP FILENAME, SEXP SIZE, SEXP SKIP);
SEXP findEOL_C(SEXP x, SEXP SKIP, SEXP HOP);
SEXP extract_POS(SEXP x, SEXP EOL, SEXP nTABtoSKIP);
SEXP extract_REF(SEXP x, SEXP EOL, SEXP nTABtoSKIP);
SEXP build_factor_loci(SEXP x, SEXP N);
SEXP summary_loci_pegas(SEXP x, SEXP LOCI);

static R_CallMethodDef Call_entries[] = {
    {"read_bin_pegas", (DL_FUNC) &read_bin_pegas, 3},
    {"findEOL_C", (DL_FUNC) &findEOL_C, 3},
    {"extract_POS", (DL_FUNC) &extract_POS, 3},
    {"extract_REF", (DL_FUNC) &extract_REF, 3},
    {"build_factor_loci", (DL_FUNC) &build_factor_loci, 2},
    {"unique_haplotype_loci", (DL_FUNC) &unique_haplotype_loci, 3},
    {"summary_loci_pegas", (DL_FUNC) &summary_loci_pegas, 2},
    {NULL, NULL, 0}
};

void R_init_pegas(DllInfo *info)
{
    R_registerRoutines(info, C_entries, Call_entries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
