#include "tabixfile.h"
#include "utilities.h"

static SEXP TABIXFILE_TAG = NULL;

static const int TBX_INIT_SIZE = 32767;

static void _tabixfile_close(SEXP ext)
{
    _TABIX_FILE *tfile = TABIXFILE(ext);
    if (NULL != tfile->tabix)
        tbx_destroy(tfile->tabix);
    tfile->tabix = NULL;
    if (NULL != tfile->iter)
        tbx_itr_destroy(tfile->iter);
    tfile->iter = NULL;
    if (NULL != tfile->fp)
        hts_close(tfile->fp);
    tfile->fp = NULL;
}

static void _tabixfile_finalizer(SEXP ext)
{
    if (NULL == R_ExternalPtrAddr(ext))
        return;
    _tabixfile_close(ext);
    _TABIX_FILE *tfile = TABIXFILE(ext);
    Free(tfile);
    R_SetExternalPtrAddr(ext, NULL);
}

SEXP tabixfile_init()
{
    TABIXFILE_TAG = install("TabixFile");
    return R_NilValue;
}

SEXP tabixfile_open(SEXP filename, SEXP indexname)
{
    if (!IS_CHARACTER(filename) || 1L != Rf_length(filename))
        Rf_error("'filename' must be character(1)");
    if (!IS_CHARACTER(indexname) || 1L != Rf_length(indexname))
        Rf_error("'indexname' must be character(1)");

    _TABIX_FILE *tfile = Calloc(1, _TABIX_FILE);
    tfile->tabix = tbx_index_load(translateChar(STRING_ELT(filename, 0)));
    if (NULL == tfile->tabix) {
        Free(tfile);
        Rf_error("failed to open tabix index file");
    }
    tfile->iter = NULL;
    tfile->fp = hts_open(translateChar(STRING_ELT(filename, 0)), "r");
    if(tfile->fp == NULL) {
        Free(tfile);
        Rf_error("failed to open tabix file '%s'",
                 translateChar(STRING_ELT(filename, 0)));
    }

    SEXP ext = PROTECT(R_MakeExternalPtr(tfile, TABIXFILE_TAG, filename));
    R_RegisterCFinalizerEx(ext, _tabixfile_finalizer, TRUE);
    UNPROTECT(1);

    return ext;
}

SEXP tabixfile_close(SEXP ext)
{
    _checkext(ext, TABIXFILE_TAG, "close");
    _tabixfile_close(ext);
    return (ext);
}

SEXP tabixfile_isopen(SEXP ext)
{
    SEXP ans = ScalarLogical(FALSE);
    if (NULL != TABIXFILE(ext)) {
        _checkext(ext, TABIXFILE_TAG, "isOpen");
        if (TABIXFILE(ext)->tabix)
            ans = ScalarLogical(TRUE);
    }
    return ans;
}

SEXP index_tabix(SEXP filename, SEXP format, SEXP seq, SEXP begin, SEXP end,
                 SEXP skip, SEXP comment, SEXP zeroBased)
{
    tbx_conf_t conf = tbx_conf_gff;

    if (!IS_CHARACTER(filename) || 1L != Rf_length(filename))
        Rf_error("'filename' must be character(1)");

    const char *fname = translateChar(STRING_ELT(filename, 0));

    if (1L == Rf_length(format)) {
        const char *txt = CHAR(STRING_ELT(format, 0));
        if (strcmp(txt, "gff") == 0)
            conf = tbx_conf_gff;
        else if (strcmp(txt, "bed") == 0)
            conf = tbx_conf_bed;
        else if (strcmp(txt, "sam") == 0)
            conf = tbx_conf_sam;
        else if (strcmp(txt, "vcf") == 0 || strcmp(txt, "vcf4") == 0)
            conf = tbx_conf_vcf;
        else if (strcmp(txt, "psltbl") == 0)
            conf = tbx_conf_psltbl;
        else
            Rf_error("format '%s' unrecognized", txt);
    } else {
        if (!IS_INTEGER(seq) || 1L != Rf_length(seq))
            Rf_error("'seq' must be integer(1)");
        conf.sc = INTEGER(seq)[0];
        if (!IS_INTEGER(begin) || 1L != Rf_length(begin))
            Rf_error("'begin' must be integer(1)");
        conf.bc = INTEGER(begin)[0];
        if (!IS_INTEGER(end) || 1L != Rf_length(end))
            Rf_error("'end' must be integer(1)");
        conf.ec = INTEGER(end)[0];
    }

    if (IS_INTEGER(skip) && 1L == Rf_length(skip))
        conf.line_skip = INTEGER(skip)[0];
    if (IS_CHARACTER(comment) && 1L == Rf_length(comment))
        conf.meta_char = CHAR(STRING_ELT(comment, 0))[0];
    if (IS_LOGICAL(zeroBased) && 1L == Rf_length(zeroBased) &&
        TRUE == LOGICAL(zeroBased)[0])
        conf.preset |= TBX_UCSC;

    if (1 != bgzf_is_bgzf(fname))
        Rf_error("file does not appear to be bgzip'd");
    if (-1L == tbx_index_build(fname, 0, &conf))
        Rf_error("index build failed");

    return filename;
}

/* returns NULL as stop condition */
const char *_tabix_read(htsFile *fp, tbx_t *t, hts_itr_t *iter)
{
    kstring_t str = {0,0,0};
    if(tbx_itr_next(fp, t, iter, &str) < 0) {
        if(str.s != NULL)
            free(str.s);
        return NULL;
    }
        /* Rf_error("error reading from tabix file '%s'", fp->fn); */
    return ks_release(&str);
}

/* Strategy change: tabix files don't respond well to bgzf_seek on
 * fp.bgzf because it's only a BGZF* some of the time. hts_open on
 * CRAM, text, sam, and vcf files all use something other than the
 * bgzf member of the htsFile::fp union member. Therefore, (1) open
 * the file in question again, (2) read the header lines, and (3)
 * close--never affecting the original htsFile parameter */
SEXP _header_lines(htsFile *fp_orig, const tbx_conf_t * conf)
{
    const int GROW_BY = 100;
    SEXP lns;
    int i_lns = 0, pidx;

    htsFile *fp = hts_open(fp_orig->fn, "r");
    if (fp == NULL)
        Rf_error("failed to open tabix file to read header lines");

    PROTECT_WITH_INDEX(lns = NEW_CHARACTER(0), &pidx);

    kstring_t str = {0, 0, 0};
    while (hts_getline(fp, '\n', &str) >= 0) {
        if (str.l != 0 && str.s[0] != conf->meta_char)
            break;
        if (0 == (i_lns % GROW_BY)) {
            lns = Rf_lengthgets(lns, Rf_length(lns) + GROW_BY);
            REPROTECT(lns, pidx);
        }
        SET_STRING_ELT(lns, i_lns++, mkCharLen(str.s, str.l));
    }

    free(str.s);
    lns = Rf_lengthgets(lns, i_lns);
    UNPROTECT(1);
    hts_close(fp);

    return lns;
}

SEXP header_tabix(SEXP ext)
{
    _checkext(ext, TABIXFILE_TAG, "scanTabix");
    tbx_t *tabix = TABIXFILE(ext)->tabix;
    htsFile *hts = TABIXFILE(ext)->fp;
    if(tabix == NULL)
        Rf_error("internal: tabix NULL in 'headerTabix' for file '%s'",
            hts->fn);

    SEXP result = PROTECT(NEW_LIST(5)), tmp, nms;
    nms = NEW_CHARACTER(Rf_length(result));
    Rf_namesgets(result, nms);
    SET_STRING_ELT(nms, 0, mkChar("seqnames"));
    SET_STRING_ELT(nms, 1, mkChar("indexColumns"));
    SET_STRING_ELT(nms, 2, mkChar("skip"));
    SET_STRING_ELT(nms, 3, mkChar("comment"));
    SET_STRING_ELT(nms, 4, mkChar("header"));

    /* seqnames */
    int n;
    const char **seqnames = tbx_seqnames(tabix, &n);
    if (n < 0)
        Rf_error("'seqnamesTabix' found < 0 (!) seqnames");
    tmp = NEW_CHARACTER(n);
    SET_VECTOR_ELT(result, 0, tmp);
    for (int i = 0; i < n; ++i)
        SET_STRING_ELT(tmp, i, mkChar(seqnames[i]));
    free(seqnames);

    const tbx_conf_t *conf = &(tabix->conf);

    /* indexColumns */
    tmp = NEW_INTEGER(3);
    SET_VECTOR_ELT(result, 1, tmp);
    INTEGER(tmp)[0] = conf->sc;
    INTEGER(tmp)[1] = conf->bc;
    INTEGER(tmp)[2] = conf->ec;
    nms = NEW_CHARACTER(3);
    Rf_namesgets(tmp, nms);
    SET_STRING_ELT(nms, 0, mkChar("seq"));
    SET_STRING_ELT(nms, 1, mkChar("start"));
    SET_STRING_ELT(nms, 2, mkChar("end"));

    /* skip */
    SET_VECTOR_ELT(result, 2, ScalarInteger(conf->line_skip));

    /* comment */
    char comment[2];
    comment[0] = (char) conf->meta_char;
    comment[1] = '\0';
    SET_VECTOR_ELT(result, 3, ScalarString(mkChar(comment)));

    /* header lines */
    SET_VECTOR_ELT(result, 4, _header_lines(hts, conf));

    UNPROTECT(1);
    return result;
}

SEXP tabix_as_character(htsFile *fp, tbx_t *tabix, hts_itr_t *iter,
                        const int size, SEXP state, SEXP rownames)
{
    const double SCALE = 1.6;
    const tbx_conf_t *conf = &(tabix->conf);

    size_t buflen = 4096;
    char *buf = Calloc(buflen, char);

    const char *line;

    int irec = 0, nrec, pidx;
    nrec = NA_INTEGER == size ? TBX_INIT_SIZE : size;
    SEXP record = NEW_CHARACTER(nrec);
    PROTECT_WITH_INDEX(record, &pidx);

    if (R_NilValue != rownames)
        Rf_error("[internal] expected 'NULL' rownames in tabix_as_character");
    if (R_NilValue != state)
        Rf_error("[internal] expected 'NULL' state in tabix_as_character");

    while (NULL != (line = _tabix_read(fp, tabix, iter)))
    {
        if (conf->meta_char == *line) {
            free((char*) line);
            line = NULL;
            continue;
        }

        if (irec == nrec) {
            nrec = nrec * SCALE;
            record = Rf_lengthgets(record, nrec);
            REPROTECT(record, pidx);
        }

        size_t linelen = strlen(line);
        if (linelen + 1 > buflen) {
            Free(buf);
            buflen = 2 * linelen;
            buf = Calloc(buflen, char);
        }
        memcpy(buf, line, linelen);
        buf[linelen] = '\0';

        SET_STRING_ELT(record, irec, mkChar(buf));

        irec += 1;
        free((char*) line);
        line = NULL;

        if (NA_INTEGER != size && irec == nrec)
            break;
    }
    if(line != NULL) {
        free((char*) line);
        line = NULL;
    }

    Free(buf);
    record = Rf_lengthgets(record, irec);
    UNPROTECT(1);
    return record;
}

SEXP tabix_count(htsFile *fp, tbx_t *tabix, hts_itr_t *iter, const int size, 
                 SEXP state, SEXP rownames)
{
    /* suppress warning about unused parameter; want consistent
     * signature for SCAN_FUN */
    (void) size;
    const tbx_conf_t *conf = &(tabix->conf);
    const char *line;
    int irec = 0;

    if (R_NilValue != rownames)
        Rf_error("[internal] expected 'NULL' rownames in tabix_count");
    if (R_NilValue != state)
        Rf_error("[internal] expected 'NULL' state in tabix_count");

    while (NULL != (line = _tabix_read(fp, tabix, iter)))
    {
        if (conf->meta_char == *line) {
            free((char*) line);
            line = NULL;
            continue;
        }
        irec += 1;
        free((char*) line);
        line = NULL;
    }
    if(line != NULL) {
        free((char*) line);
        line = NULL;
    }

    return ScalarInteger(irec);
}

SEXP scan_tabix(SEXP ext, SEXP space, SEXP yield, SEXP fun,
                SEXP state, SEXP rownames)
{
    _checkparams(space, R_NilValue, R_NilValue);
    if (!IS_INTEGER(yield) || 1L != Rf_length(yield))
        Rf_error("'yieldSize' must be integer(1)");
    _checkext(ext, TABIXFILE_TAG, "scanTabix");

    tbx_t *tabix = TABIXFILE(ext)->tabix;
    htsFile *fp = TABIXFILE(ext)->fp;
    if(tabix == NULL || fp == NULL) {
        Rf_error("internal: tabix or htsFile missing; not properly opened?");
    }
    SCAN_FUN *scan = (SCAN_FUN *) R_ExternalPtrAddr(fun);

    SEXP spc = VECTOR_ELT(space, 0);
    const int nspc = Rf_length(spc);
    SEXP result, elt;

    PROTECT(result = NEW_LIST(nspc == 0 ? 1 : nspc));
    if (0 == nspc) {
        hts_itr_t *iter = TABIXFILE(ext)->iter;
        if (NULL == iter) {
            /* HTS_IDX_START doesn't work because there aren't aligned
             * reads(?); HTS_IDX_REST signals to iterate from the
             * current offset over the rest of the file */
            iter = TABIXFILE(ext)->iter = tbx_itr_queryi(tabix,
                                                         HTS_IDX_REST, -1, -1);
            if(iter == NULL)
                Rf_error("problem creating iterator for file '%s'", fp->fn);
        }
        elt = scan(fp, tabix, iter, INTEGER(yield)[0],
                   state, rownames);
        SET_VECTOR_ELT(result, 0, elt);
    } else {
        const int
            *start = INTEGER(VECTOR_ELT(space, 1)),
            *end = INTEGER(VECTOR_ELT(space, 2));

        for (int ispc = 0; ispc < nspc; ++ispc) {
            int ibeg, iend, tid;
            hts_itr_t *iter;
            const char *tid_name;

            ibeg = start[ispc] == 0 ? 0 : start[ispc] - 1;
            iend = end[ispc];
            tid_name = CHAR(STRING_ELT(spc, ispc));
            /* char regbuf[256]; */
            /* sprintf(regbuf, "%s:%d-%d", tid_name, ibeg, iend); */
            tid = tbx_name2id(tabix, tid_name);
            if (0 > tid)
                Rf_error("'%s' not present in tabix index", tid_name);
            iter = tbx_itr_queryi(tabix, tid, ibeg, iend);
            /* iter = tbx_itr_querys(tabix, regbuf); */

            elt = scan(fp, tabix, iter, NA_INTEGER,
                       state, rownames);
            SET_VECTOR_ELT(result, ispc, elt);

            tbx_itr_destroy(iter);
        }
    }

    UNPROTECT(1);
    return result;
}
