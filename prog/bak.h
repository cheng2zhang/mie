/* compute the reference entropy */
__inline static double potts_ent2ref(potts_t *p, double beta)
{
  int s, t, st, q = p->q, i, j, n = p->n;
  double **lnm; /* transition matrices, saving logarithms of elements */
  double *lnz; /* partition function */
  double lnp, eij;

  /* allocate spaces for transition matrices */
  xnew(lnm, n);
  for ( i = 0; i < n; i++ ) {
    xnew(lnm[i], q * q);
  }
  xnew(lnz, n);

  for ( s = 0; s < q; s++ ) {
    for ( t = 0; t < q; t++ ) {
      st = s * q + t;
      lnm[0][st] = beta * (s == t); /* exp(-beta*(-delta(s,t))) */
      lnm[1][st] = lnm[0][st];
    }
  }

  /* using matrix multiplication
   * to compute higher-order matrices */
  for ( i = 1; i < n - 1; i++ ) {
    /* lnm[i+1] = lnm[i] * lnm[0]; */
    lnmatmul(lnm[i+1], lnm[i], lnm[0], q);
    fprintf(stderr, "matrix %d: %g, %g\n", i+1, exp(lnm[i+1][0]), exp(lnm[i+1][1]));
  }

  /* compute the partition function */
  for ( i = 0; i < n; i++ ) {
    lnz[i] = LN0;
    for ( s = 0; s < q; s++ ) {
      for ( t = 0; t < q; t++ ) {
        st = s * q + t;
        lnz[i] = lnadd(lnz[i], lnm[i][st]);
      }
    }
    if ( i > 0 ) {
      fprintf(stderr, "i %d: %g, %g\n", i, lnz[i], log(q) + i*log(q+exp(beta)-1));
    }
  }

  /* compute the information */
  p->ent2r = 0;
  for ( i = 0; i < n; i++ ) {
    for ( j = i + 1; j < n; j++ ) {
      eij = 0;
      for ( s = 0; s < q; s++ ) {
        for ( t = 0; t < q; t++ ) {
          st = s * q + t;
          lnp = lnm[j-i][st] - lnz[j-i];
          eij += -exp(lnp) * lnp;
        }
      }
      p->ent2r += eij - 2*log(q);
    }
  }
  p->ent1r = n * log(q);
  p->ent2r += p->ent1r;

  /* compute the pairwise information */
  for ( i = 0; i < n; i++ ) free(lnm[i]);
  free(lnm);
  free(lnz);

  return p->ent2r;
}
