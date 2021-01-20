#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(zooem)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zoolayerem)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zoofb)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zoolayerfb)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zoolik)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zoolayerlik)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zooviterbi)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zoolayerviterbi)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zoosim)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void*, void*);

static const R_FortranMethodDef FortranEntries[] = {
    {"zooem",      (DL_FUNC) &F77_NAME(zooem),      17},
    {"zoolayerem", (DL_FUNC) &F77_NAME(zoolayerem),      17},
    {"zoofb",      (DL_FUNC) &F77_NAME(zoofb),      10},
    {"zoolayerfb", (DL_FUNC) &F77_NAME(zoolayerfb),      10},
    {"zoolik",     (DL_FUNC) &F77_NAME(zoolik),      9},
    {"zoolayerlik",(DL_FUNC) &F77_NAME(zoolayerlik),      9},
    {"zooviterbi", (DL_FUNC) &F77_NAME(zooviterbi),  9},
    {"zoolayerviterbi", (DL_FUNC) &F77_NAME(zoolayerviterbi),  9},
    {"zoosim",     (DL_FUNC) &F77_NAME(zoosim),     13},
    {NULL, NULL, 0}
};

void R_init_RZooRoH(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

