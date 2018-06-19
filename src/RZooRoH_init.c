#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(zooem)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zoofb)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zoolik)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zooviterbi)(void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"zooem",      (DL_FUNC) &F77_NAME(zooem),      17},
    {"zoofb",      (DL_FUNC) &F77_NAME(zoofb),      10},
    {"zoolik",     (DL_FUNC) &F77_NAME(zoolik),      9},
    {"zooviterbi", (DL_FUNC) &F77_NAME(zooviterbi),  9},
    {NULL, NULL, 0}
};

void R_init_RZooRoH(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

