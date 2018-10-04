#define IS2_LB 4
#define L  (1 << IS2_LB)
#include "is2.h"
#include "ave.h"
#include "argopt.h"
#include <time.h>

double tp = 2.2;
int nsys = 100;
long nsteps = 10000;
int method = 0; /* 0: Metropolis, 1: Wolff */
long nequil = 10*IS2_N;
int npart = 2;
int npart2 = 4;
long nstrep = 100;
char *fnlog = "is2.log";

static void doargs(int argc, char **argv)
{
  argopt_t *ao;

  ao = argopt_open(0);
  ao->desc = "Entropy estimation along magnetization";
  argopt_add(ao, "-T", "%lf", &tp, "temperature");
  argopt_add(ao, "-M", "%d", &nsys, "number of replica systems");
  argopt_add(ao, "-t", "%ld", &nsteps, "number of steps");
  argopt_add(ao, "-r", "%ld", &nstrep, "number of steps to report");
  argopt_add(ao, "-P", "%d", &npart, "number of partitions for the block method");
  argopt_add(ao, "-Q", "%d", &npart2, "number of partitions for the block method (2)");
  argopt_add(ao, "-m", "%d", &method, "sampling method");
  argopt_add(ao, "-E", "%ld", &nequil, "number of equilibration");
  argopt_add(ao, "-g", NULL, &fnlog, "name of the log file");
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);
  argopt_dump(ao);
  argopt_close(ao);
}

typedef struct {
  is2_t *is;
  int q;
  int trjn, trji, *trj;
  long tot, *cnt; /* histogram data */
  double ent, entc, entcb, ents, entsb;
  double mest;
} walker_t;

static walker_t *walker_open(long nsteps)
{
  walker_t *w;
  int s, q;

  xnew(w, 1);
  w->is = is2_open(L);
  w->q = q = 2*w->is->n + 1;
  w->tot = 0;
  xnew(w->cnt, q);
  for ( s = 0; s < q; s++ )
    w->cnt[s] = 0;

  /* initialize the trajectory data */
  w->trjn = nsteps;
  w->trji = 0;
  xnew(w->trj, w->trjn);
  return w;
}

/* count occurrences from frame start to frame end */
static void walker_count(walker_t *w, int start, int end)
{
  int i, s;

  for ( s = 0; s < w->q; s++ )
    w->cnt[s] = 0;
  for ( i = start; i < end; i++ ) {
    w->cnt[ w->trj[i] ] += 1;
  }
}


static double walker_ent_seg(walker_t *w, int start, int end)
{
  double ent = 0, tot = (double) (end - start), pr;
  int s;

  walker_count(w, start, end);
  for ( s = 0; s < w->q; s++ ) {
    long c = w->cnt[s];
    if ( c <= 0 ) continue;
    pr = 1.0*c/tot;
    ent += -pr*log(pr);
  }
  return ent;
}

static double walker_entropy(walker_t *w, int npart, int npart2, int verbose)
{
  int trjn = w->trji, ip, blksz, blksz2;
  double entp = 0, entp2 = 0, entt, entc2;

  /* entropy estimate from the entire trajectory */
  entt = walker_ent_seg(w, 0, trjn);

  /* entropy estimated from trajectory block(s) */
  blksz = trjn / npart;
  for ( ip = 0; ip < npart; ip++ ) {
    entp += walker_ent_seg(w, ip * blksz, (ip + 1) * blksz);
  }
  entp /= npart;

  blksz2 = trjn / npart2;
  for ( ip = 0; ip < npart2; ip++ ) {
    entp2 += walker_ent_seg(w, ip * blksz2, (ip + 1) * blksz2);
  }
  entp2 /= npart2;

  /* since we have St = S - a/t, Sp = S - a/blksz
   * S = (St*t - Sp*blksz)/(t-blksz); */
  w->ent = entt;
  w->entc = (entt*trjn - entp*blksz)/(trjn - blksz);
  entc2 = (entt*trjn - entp2*blksz2)/(trjn - blksz2);
  w->entcb = (w->entc*blksz - entc2*blksz2)/(blksz - blksz2);

  /* fitting to log(1+m/2/t) */
  {
    double ds = entt - entp, xp = exp(ds);
    double ds2 = entt - entp2, xp2 = exp(ds2);
    double xpc, xpc2, xpcb;

    {
      if ( xp >= .99*trjn/blksz ) xp = .99*trjn/blksz;
      xpc = (trjn - xp*blksz)/(trjn - blksz);
      w->ents = entt - log(xpc);
      w->mest = 2*trjn*(xp - 1)/(1.*trjn/blksz - xp);

      if ( xp2 >= .99*trjn/blksz2 ) xp2 = .99*trjn/blksz2;
      xpc2 = (trjn - xp2*blksz2)/(trjn - blksz2);

      xpcb = (xpc*blksz - xpc2*blksz2)/(blksz - blksz2);
      if ( xpcb < 0 ) xpcb = xpc;
      w->entsb = entt - log(xpcb);
    }
  }

  return w->ent;
}

static void walker_close(walker_t *w)
{
  is2_close(w->is);
  free(w->cnt);
  free(w->trj);
  free(w);
}

static void work(void)
{
  walker_t **w;
  av_t avent[1], aventc[1], aventcb[1], avents[1], aventsb[1], avmest[1];
  double ent, entc, entcb, ents, entsb, mest;
  double var, varc, varcb, vars, varsb, mvar;
  int i;
  long t;
  FILE *fplog;

  int id, h;
  double beta = 1.0/tp, padd;
  double entref = 0;

  if ( (fplog = fopen(fnlog, "w")) == NULL ) {
    fprintf(stderr, "cannot open log file %s\n", fnlog);
    fplog = stderr;
  }

  padd = 1 - exp(-2*beta);
  xnew(w, nsys);
  for ( i = 0; i < nsys; i++ ) {
    w[i] = walker_open(nsteps);
    is2_setuproba(beta, w[i]->is->uproba);
    /* equilibration */
    for ( t = 1; t <= nequil; t++ ) {
      is2_wolff(w[i]->is, padd);
    }
  }

  /* production run */
  for ( t = 1; t <= nsteps; t++ ) {
    for ( i = 0; i < nsys; i++ ) {
      is2_t *is = w[i]->is;
  
      if ( method == 0 ) { /* Metropolis algorithm */
        //id = is2_pick(is, &h);
        IS2_PICK(is, id, h);
        if ( h <= 0 || mtrand() <= is->uproba[h] ) {
          //is2_flip(is, id, h);
          IS2_FLIP(is, id, h);
        }
      } else { /* Wolff cluster algorithm */
        is2_wolff(is, padd);
      }

      w[i]->trj[ w[i]->trji++ ] = is->M + IS2_N;
    }

    if ( t % nstrep == 0 ) {
      av_clear(avent);
      av_clear(aventc);
      av_clear(aventcb);
      av_clear(avents);
      av_clear(aventsb);
      av_clear(avmest);
      for ( i = 0; i < nsys; i++ ) {
        walker_entropy(w[i], npart, npart2, i <= 0);
        av_add(avent,  w[i]->ent);
        av_add(aventc, w[i]->entc);
        av_add(aventcb, w[i]->entcb);
        av_add(avents, w[i]->ents);
        av_add(aventsb, w[i]->entsb);
        av_add(avmest, w[i]->mest);
      }
      ent  = av_getave(avent,  &var);
      entc = av_getave(aventc, &varc);
      entcb = av_getave(aventcb, &varcb);
      ents = av_getave(avents, &vars);
      entsb = av_getave(aventsb, &varsb);
      mest = av_getave(avmest, &mvar);
      printf("%9ld: entropy %8.4f(%g), %8.4f, %8.4f; %8.4f, %8.4f(%8.4f), m %6.0f(%5.0f)\n",
          t, ent, var, entc, entcb, ents, entsb, entref, mest, sqrt(mvar));
      fprintf(fplog, "%ld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", t,
          ent, sqrt(var), entref,
          entc, sqrt(varc), entcb, sqrt(varcb),
          ents, sqrt(vars), entsb, sqrt(varsb),
          mest, sqrt(mest));
      fflush(fplog);
    }
  }

  for ( i = 0; i < nsys; i++ ) {
    walker_close(w[i]);
  }
  free(w);
  fclose(fplog);
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  work();
  return 0;
}
