#include "../mtrand.h"
#include "../ave.h"
#include "../argopt.h"
#include <time.h>

int blksz = 1000;
char *fnseq = "../../data/words/shakespeare/shakespeare.seq";
char *fnlog = "a.log";

static void doargs(int argc, char **argv)
{
  argopt_t *ao;

  ao = argopt_open(0);
  ao->desc = "Word walk";
  argopt_add(ao, "-i", NULL,  &fnseq,     "sequence file");
  argopt_add(ao, "-b", "%d",  &blksz,     "block size");
  argopt_add(ao, "-o", "%s",  &fnlog,     "output log file");
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);
  argopt_dump(ao);
  argopt_close(ao);
}

typedef struct {
  int nw;
  int nseq;
  int **seq;
  int *len;
} seq_t;


static seq_t *loadseq(const char *fn)
{
  FILE *fp;
  char buf[1024];
  int i, j, id, len;
  seq_t *seq;

  xnew(seq, 1);

  /* load the sequence */
  if ((fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fnseq);
    return NULL;
  }
  fgets(buf, sizeof buf, fp);
  sscanf(buf + 1, "%d %d", &seq->nw, &seq->nseq);
  xnew(seq->seq, seq->nseq);
  xnew(seq->len, seq->nseq);
  for ( i = 0; i < seq->nseq; i++ ) {
    fgets(buf, sizeof buf, fp);
    sscanf(buf + 1, "%d %d", &id, &len);
    if ( id != i ) {
      fprintf(stderr, "sequence identity corruption\n");
      exit(1);
    }
    seq->len[i] = len;
    xnew(seq->seq[i], len);
    for ( j = 0; j < len; j++ ) {
      fgets(buf, sizeof buf, fp);
      sscanf(buf, "%d", &seq->seq[i][j]);
    }
    fgets(buf, sizeof buf, fp); /* a blank line */
  }
  fclose(fp);
  fprintf(stderr, "got %d sequences of %d words\n", seq->nseq, seq->nw);

  return seq;
}

/* entropy by counting */
static double get_seq_entropy(const int *seq, int len, int *cnt, int ntok)
{
  double ent = 0, x = 0.0;
  int i, s;

  for ( s = 0; s < ntok; s++ ) {
    cnt[s] = 0;
  }
  for ( i = 0; i < len; i++ ) {
    s = seq[i];
    cnt[s] += 1;
  }
  for ( s = 0; s < ntok; s++ ) {
    if ( cnt[s] <= 0 ) continue;
    x = 1.0*cnt[s]/len;
    ent -= x*log(x);
  }
  return ent;
}

static double get_seq_entropy_extra(const int *seq, int len,
    int *cnt, int ntok, double *entlin, double *entexp)
{
  double ent = 0, ent1 = 0;
  int blksz;

  ent = get_seq_entropy(seq, len, cnt, ntok);

  blksz = len/2;
  ent1  = get_seq_entropy(seq, blksz, cnt, ntok);
  ent1 += get_seq_entropy(seq + blksz, len - blksz, cnt, ntok);
  ent1 /= 2;

  *entlin = 2*ent - ent1;
  *entexp = ent - log(2 - exp(ent - ent1));
  return ent;
}



/*
static double get_block_entropy1(const int *seq, int len, int ntok,
    int blksz, int *cnt)
{
  double ent = 0, ent1 = 0;
  int ib, nb;

  nb = len / blksz;
  for ( ib = 0; ib < nb; ib++ ) {
    ent1 = get_seq_entropy(seq + ib*blksz, blksz, cnt, ntok);
    ent += ent1;
    //printf("%d %g %g\n", ib, ent1, ent/(ib+1));
  }
  ent = ent / nb;
  return ent;
}


static double get_block_entropy2(const int *seq, int len, int ntok,
    int blksz, int *cnt, int del)
{
  double ent = 0, ent1 = 0;
  int i0, tot = 0;

  for ( i0 = 0; i0 + blksz <= len; i0 += del ) {
    ent1 = get_seq_entropy(seq + i0, blksz, cnt, ntok);
    ent += ent1;
    tot += 1;
  }
  ent = ent / tot;
  return ent;
}
*/
int work(void)
{
  seq_t *seq;
  int *cnt = NULL, i;
  double ent, entlin, entexp, entref;
  FILE *fp;

  seq = loadseq(fnseq);
  xnew(cnt, seq->nw); /* counter */

  if ((fp = fopen(fnlog, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fnlog);
    return -1;
  }

  for ( i = 0; i < seq->nseq; i++ ) {
    ent = get_seq_entropy_extra(seq->seq[i], seq->len[i], cnt, seq->nw, &entlin, &entexp);
    fprintf(fp, "%d %g %g %g\n", seq->len[i], ent, entlin, entexp);
    printf("%d %g %g %g\n", seq->len[i], ent, entlin, entexp);
  }

  fclose(fp);
  fprintf(stderr, "saved result to %s\n", fnlog);
  free(cnt);
  free(seq);
  return 0;
}


int main(int argc, char **argv)
{
  doargs(argc, argv);
  work();
  return 0;
}
