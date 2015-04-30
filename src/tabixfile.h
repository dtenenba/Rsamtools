#ifndef _TABIXFILE_H
#define _TABIXFILE_H

#include <Rdefines.h>
#include "tbx.h"
#include "bgzf.h"
#include "kstring.h"

typedef struct {
    tbx_t *tabix;
    hts_itr_t *iter;
    htsFile *fp;
} _TABIX_FILE;

#define TABIXFILE(b) ((_TABIX_FILE *) R_ExternalPtrAddr(b))

typedef SEXP SCAN_FUN(htsFile *fp, tbx_t *tabix, hts_itr_t *iter, const int size,
                      SEXP state, SEXP rownames);

SCAN_FUN tabix_as_character;
SCAN_FUN tabix_count;

SEXP tabixfile_init();
SEXP tabixfile_open(SEXP filename, SEXP indexname);
SEXP tabixfile_close(SEXP ext);
SEXP tabixfile_isopen(SEXP ext);

SEXP index_tabix(SEXP filename, SEXP format,
                 SEXP seq, SEXP begin, SEXP end,
                 SEXP skip, SEXP comment, SEXP zeroBased);
SEXP header_tabix(SEXP ext);
SEXP scan_tabix(SEXP ext, SEXP space, SEXP yield, SEXP fun, 
                SEXP state, SEXP rownames);

#endif
