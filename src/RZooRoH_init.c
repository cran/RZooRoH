#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(zoolayerfb)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zoolayerlik)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zoosumlayerlik)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zoosumlayerlik2)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zoosumlayerfb)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void*);
extern void F77_NAME(zoosumlayerfb2)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void*);
extern void F77_NAME(zoolayerviterbi)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zoosumlayerviterbi)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zoosumlayerviterbi2)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(freqem1)(void *, void *, void *);
extern void F77_NAME(freqem2)(void *, void *, void *);
extern void F77_NAME(freqem3)(void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"zoolayerfb", (DL_FUNC) &F77_NAME(zoolayerfb),      10},
    {"zoolayerlik",(DL_FUNC) &F77_NAME(zoolayerlik),      9},
    {"zoosumlayerlik",(DL_FUNC) &F77_NAME(zoosumlayerlik),      9},
    {"zoosumlayerlik2",(DL_FUNC) &F77_NAME(zoosumlayerlik2),      9},
    {"zoosumlayerfb",(DL_FUNC) &F77_NAME(zoosumlayerfb),      10},
    {"zoosumlayerfb2",(DL_FUNC) &F77_NAME(zoosumlayerfb2),      10},
    {"zoolayerviterbi", (DL_FUNC) &F77_NAME(zoolayerviterbi),  9},
    {"zoosumlayerviterbi", (DL_FUNC) &F77_NAME(zoosumlayerviterbi),  9},
    {"zoosumlayerviterbi2", (DL_FUNC) &F77_NAME(zoosumlayerviterbi2),  9},
    {"freqem1",     (DL_FUNC) &F77_NAME(freqem1),     3},
    {"freqem2",     (DL_FUNC) &F77_NAME(freqem2),     3},
    {"freqem3",     (DL_FUNC) &F77_NAME(freqem3),     3},
    {NULL, NULL, 0}
};

void R_init_RZooRoH(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

