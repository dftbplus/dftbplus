#ifndef blacsutils_h
#define blacsutils_h

#define    DTYPE_             0                   /* Descriptor Type */
#define    CTXT_              1                     /* BLACS context */
#define    M_                 2             /* Global Number of Rows */
#define    N_                 3          /* Global Number of Columns */
#define    MB_                4                 /* Row Blocking Size */
#define    NB_                5              /* Column Blocking Size */
#define    RSRC_              6            /* Starting Processor Row */
#define    CSRC_              7         /* Starting Processor Column */
#define    LLD_               8           /* Local Leading Dimension */
#define    DLEN_              9                 /* Descriptor Length */

extern void blacs_pinfo_(int *mypnum, int *nprocs);
extern void blacs_get_(int *, int *, int *);
extern void blacs_gridinit_(
    int *, char *, int *,
    int *);  // BLACS_GRIDINIT( ICONTXT, ORDER, NPROW, NPCOL )

// https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/scalapack-routines/scalapack-redistribution-copy-routines/p-gemr2d.html
// http://www.netlib.org/scalapack/slug/node164.html
// http://www.netlib.org/scalapack/slug/node168.html
extern void pdgemr2d_(int *m, int *n, double *a, int *ia, int *ja, int *desca,
                      double *b, int *ib, int *jb, int *descb, int *ictxt);
extern void descinit_(
    int *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *);  // DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT, LLD, INFO

extern void blacs_gridinfo_(
    int *, int *, int *, int *,
    int *);  // blacs_gridinfo (icontxt, nprow, npcol, myrow, mycol);
extern int numroc_(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);

int get_system_context(int blacs_context) {
  int what = 10, val;
  blacs_get_(&blacs_context, &what,
             &val);  // get system context of the blacs_context
  return val;
}

int get_default_system_context() {
  int system_context, zero = 0;
  blacs_get_(&zero, &zero, &system_context);  // get default system context
  return system_context;
}

int make_blacs_context(int init_context, int MP, int NP) {
  char order = 'R';
  blacs_gridinit_(&init_context, &order, &MP, &NP);

  return init_context;
}

void blacs_desc_init(int M, int N, int blacs_context, int *desc) {
  if (blacs_context == -1) {
    desc[CTXT_] = -1;
    return;  // stub for non-participating process
  }

  int nprow, npcol, myrow, mycol;
  blacs_gridinfo_(&blacs_context, &nprow, &npcol, &myrow, &mycol);

  int MB = 1, NB = 1, rsrc = 0, csrc = 0;
  int info;

  int lld = numroc_(&M, &MB, &myrow, &rsrc, &nprow);
  if (lld < 1) {
    lld = 1;
  }

  descinit_(desc, &M, &N, &MB, &NB, &rsrc, &csrc, &blacs_context, &lld,
            &info);  // DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT, LLD, INFO
}

int locrow(int *desc) {
  int nprow, npcol, myrow, mycol;
  blacs_gridinfo_(&desc[CTXT_], &nprow, &npcol, &myrow, &mycol);

  return numroc_(&desc[M_], &desc[MB_], &myrow, &desc[RSRC_], &nprow);
}

int loccol(int *desc) {
  int nprow, npcol, myrow, mycol;
  blacs_gridinfo_(&desc[CTXT_], &nprow, &npcol, &myrow, &mycol);

  return numroc_(&desc[N_], &desc[NB_], &mycol, &desc[CSRC_], &npcol);
}

#endif  // blacsutils_h
