#include "samtools/khash.h"
#include "samtools/sam.h"
#include "bam_data.h"
#include "scan_bam_data.h"
#include "io_sam.h"
#include "bamfile.h"
#include "encode.h"
#include "utilities.h"
#include "IRanges_interface.h"
#include "XVector_interface.h"
#include "Biostrings_interface.h"
#include "bam_mate_iter.h"

/* from samtoools/bam_sort.c */
void bam_sort_core(int is_by_qname, const char *fn, const char *prefix,
                   size_t max_mem);

#define SEQUENCE_BUFFER_ALLOCATION_ERROR 1

typedef int (*_PARSE1_FUNC) (const bam1_t *, void *);
typedef void (_FINISH1_FUNC) (BAM_DATA);

static const char *TMPL_ELT_NMS[] = {
    "qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "cigar",
    "mrnm", "mpos", "isize", "seq", "qual", "tag"
    /* "vtype", "value" */
};

static const int N_TMPL_ELTS = sizeof(TMPL_ELT_NMS) / sizeof(const char *);

/* utility */

void _check_is_bam(const char *filename)
{
    int magic_len;
    char buf[4];
    bamFile bfile = bam_open(filename, "r");

    if (bfile == 0)
        Rf_error("failed to open SAM/BAM file\n  file: '%s'", filename);

    magic_len = bam_read(bfile, buf, 4);
    bam_close(bfile);

    if (magic_len != 4 || strncmp(buf, "BAM\001", 4) != 0)
        Rf_error("'filename' is not a BAM file\n  file: %s", filename);
}

/* template */

void _bam_check_template_list(SEXP template_list)
{
    if (!IS_LIST(template_list) || LENGTH(template_list) != N_TMPL_ELTS)
        Rf_error("'template' must be list(%d)", N_TMPL_ELTS);
    SEXP names = GET_ATTR(template_list, R_NamesSymbol);
    if (!IS_CHARACTER(names) || LENGTH(names) != N_TMPL_ELTS)
        Rf_error("'names(template)' must be character(%d)", N_TMPL_ELTS);
    for (int i = 0; i < LENGTH(names); ++i)
        if (strcmp(TMPL_ELT_NMS[i], CHAR(STRING_ELT(names, i))) != 0)
            Rf_error("'template' names do not match scan_bam_template\n'");
}

static SEXP _tmpl_DNAStringSet()
{
    CharAEAE aeae = new_CharAEAE(0, 0);
    SEXP lkup = PROTECT(_get_lkup("DNAString"));
    SEXP ans = new_XRawList_from_CharAEAE("DNAStringSet", "DNAString",
                                          &aeae, lkup);
    UNPROTECT(1);
    return ans;
}

static SEXP _tmpl_BStringSet()
{
    CharAEAE aeae = new_CharAEAE(0, 0);
    return new_XRawList_from_CharAEAE("BStringSet", "BString", &aeae,
                                      R_NilValue);
}

static SEXP _tmpl_PhredQuality()
{
    SEXP xstringset, s, t, nmspc, result;
    PROTECT(xstringset = _tmpl_BStringSet());
    PROTECT(nmspc = _get_namespace("Rsamtools"));
    NEW_CALL(s, t, "PhredQuality", nmspc, 2);
    CSET_CDR(t, "x", xstringset);
    CEVAL_TO(s, nmspc, result);
    UNPROTECT(2);
    return result;
}

SEXP scan_bam_template(SEXP tag)
{
    if (R_NilValue != tag)
        if (!IS_CHARACTER(tag))
            Rf_error("'tag' must be NULL or 'character()'");
    SEXP tmpl = PROTECT(NEW_LIST(N_TMPL_ELTS));
    SET_VECTOR_ELT(tmpl, QNAME_IDX, NEW_CHARACTER(0));
    SET_VECTOR_ELT(tmpl, FLAG_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, RNAME_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, STRAND_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, POS_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, QWIDTH_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, MAPQ_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, CIGAR_IDX, NEW_CHARACTER(0));
    SET_VECTOR_ELT(tmpl, MRNM_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, MPOS_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, ISIZE_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, SEQ_IDX, _tmpl_DNAStringSet());
    SET_VECTOR_ELT(tmpl, QUAL_IDX, _tmpl_PhredQuality());
    if (R_NilValue == tag) {
        SET_VECTOR_ELT(tmpl, TAG_IDX, R_NilValue);
    } else {
        SET_VECTOR_ELT(tmpl, TAG_IDX, NEW_LIST(LENGTH(tag)));
        SET_ATTR(VECTOR_ELT(tmpl, TAG_IDX), R_NamesSymbol, tag);
    }

    SEXP names = PROTECT(NEW_CHARACTER(N_TMPL_ELTS));
    for (int i = 0; i < N_TMPL_ELTS; ++i)
        SET_STRING_ELT(names, i, mkChar(TMPL_ELT_NMS[i]));
    SET_ATTR(tmpl, R_NamesSymbol, names);
    UNPROTECT(2);
    return tmpl;
}

/* header */

SEXP _read_bam_header(SEXP ext)
{
    samfile_t *sfile = BAMFILE(ext)->file;
    bam_header_t *header = sfile->header;
    int n_elts = header->n_targets;

    SEXP ans = PROTECT(NEW_LIST(2));

    /* target length / name */
    SET_VECTOR_ELT(ans, 0, NEW_INTEGER(n_elts));
    SEXP tlen = VECTOR_ELT(ans, 0); /* target length */
    SEXP tnm = PROTECT(NEW_CHARACTER(n_elts));  /* target name */
    setAttrib(tlen, R_NamesSymbol, tnm);
    UNPROTECT(1);
    for (int j = 0; j < n_elts; ++j) {
        INTEGER(tlen)[j] = header->target_len[j];
        SET_STRING_ELT(tnm, j, mkChar(header->target_name[j]));
    }

    /* 'aux' character string */
    char *txt = (char *) R_alloc(header->l_text + 1, sizeof(char));
    strncpy(txt, header->text, header->l_text);
    txt[header->l_text] = '\0';
    SET_VECTOR_ELT(ans, 1, mkString(txt));

    SEXP nms = PROTECT(NEW_CHARACTER(2));
    SET_STRING_ELT(nms, 0, mkChar("targets"));
    SET_STRING_ELT(nms, 1, mkChar("text"));
    setAttrib(ans, R_NamesSymbol, nms);
    UNPROTECT(2);
    return ans;
}

/* scan_bam framework */

static int _scan_bam_all(BAM_DATA bd, _PARSE1_FUNC parse1,
                         _FINISH1_FUNC finish1)
{
    bam1_t *bam = bam_init1();
    int r = 0;
    BAM_FILE bfile = _bam_file_BAM_DATA(bd);

    bam_seek(bfile->file->x.bam, bfile->pos0, SEEK_SET);
    int qname_bufsize = 1000;
    char *last_qname = Calloc(qname_bufsize, char);
    int yieldSize = bd->yieldSize,
        obeyQname = bd->obeyQname,
        asMates = bd->asMates;
    int ith_yield = 0, inc_yield = 1;

    while ((r = samread(bfile->file, bam)) >= 0) {
        if (NA_INTEGER != yieldSize) {
            if (obeyQname) {
                if (0 != strcmp(last_qname, bam1_qname(bam))) {
                    inc_yield = 1;
                    if (ith_yield >= yieldSize)
                        break;
                    if (bam->core.l_qname > qname_bufsize) {
                        Free(last_qname);
                        qname_bufsize = bam->core.l_qname;
                        last_qname = Calloc(qname_bufsize, char);
                    }
                    strcpy(last_qname, bam1_qname(bam));
                } else {
                    inc_yield = 0;
                }
            }
        }

        int result = (*parse1) (bam, bd);
        if (result < 0) {   /* parse error: e.g., cigar buffer overflow */
            _grow_SCAN_BAM_DATA(bd, 0);
            return result;
        } else if (result == 0L) /* does not pass filter */
            continue;

        ith_yield += inc_yield;
        if (NA_INTEGER != yieldSize && ith_yield == yieldSize) {
            bfile->pos0 = bam_tell(bfile->file->x.bam);
            if (!obeyQname)
                break;
        }
    }

    if (NULL != finish1)
        (*finish1) (bd);
    if ((NA_INTEGER == yieldSize) || (ith_yield < yieldSize))
        /* end-of-file */
        bfile->pos0 = bam_tell(bfile->file->x.bam);

    Free(last_qname);
    return bd->iparsed;
}

static int _scan_bam_fetch(BAM_DATA bd, SEXP space, int *start, int *end,
                           _PARSE1_FUNC parse1, _FINISH1_FUNC finish1)
{
    int tid;
    BAM_FILE bfile = _bam_file_BAM_DATA(bd);
    samfile_t *sfile = bfile->file;
    bam_index_t *bindex = bfile->index;
    int n_tot = bd->iparsed,
        asMates = bd->asMates;

    for (int irange = 0; irange < LENGTH(space); ++irange) {
        const char *spc = translateChar(STRING_ELT(space, irange));
        const int starti =
            start[irange] > 0 ? start[irange] - 1 : start[irange];
        for (tid = 0; tid < sfile->header->n_targets; ++tid) {
            if (strcmp(spc, sfile->header->target_name[tid]) == 0)
                break;
        }
        if (tid == sfile->header->n_targets) {
            Rf_warning("space '%s' not in BAM header", spc);
            return -1;
        }
        if (asMates)
            bam_mate_fetch(sfile->x.bam, bindex, tid, starti, 
                           end[irange], bd, parse1);
        else
            bam_fetch(sfile->x.bam, bindex, tid, starti, 
                      end[irange], bd, parse1);
        if (NULL != finish1)
            (*finish1) (bd);
        bd->irange += 1;
    }
    return bd->iparsed - n_tot;
}

static int _do_scan_bam(BAM_DATA bd, SEXP space, _PARSE1_FUNC parse1,
                        _FINISH1_FUNC finish1)
{
    int status;

    if (R_NilValue == space)
        /* everything */
        status = _scan_bam_all(bd, parse1, finish1);
    else {                      /* fetch */
        BAM_FILE bfile = _bam_file_BAM_DATA(bd);
        if (NULL == bfile->index)
            Rf_error("valid 'index' file required");
        status = _scan_bam_fetch(bd, VECTOR_ELT(space, 0),
                                 INTEGER(VECTOR_ELT(space, 1)),
                                 INTEGER(VECTOR_ELT(space, 2)),
                                 parse1, finish1);
    }

    return status;
}

/* scan_bam */

static int _filter_and_parse1(const bam1_t *bam, void *data)
{
    return _filter_and_parse1_BAM_DATA(bam, (BAM_DATA) data);
}

SEXP _scan_bam_result_init(SEXP template_list, SEXP names, SEXP space)
{
    const int nrange =
        R_NilValue == space ? 1 : Rf_length(VECTOR_ELT(space, 0));
    int i;

    SEXP result = PROTECT(NEW_LIST(nrange));
    /* result:
       range1: tmpl1, tmpl2...
       range2: tmpl1, tmpl2...
       ...
    */
    for (int irange = 0; irange < nrange; ++irange) {
        SEXP tag = VECTOR_ELT(template_list, TAG_IDX);
        SEXP tmpl;
        if (R_NilValue == tag)
            tmpl = PROTECT(scan_bam_template(R_NilValue));
        else {
            SEXP nms = getAttrib(tag, R_NamesSymbol);
            tmpl = PROTECT(scan_bam_template(nms));
        }
        for (i = 0; i < LENGTH(names); ++i) {
            if (TAG_IDX == i)
                continue;
            else if (R_NilValue == VECTOR_ELT(template_list, i))
                SET_VECTOR_ELT(tmpl, i, R_NilValue);
        }
        SET_VECTOR_ELT(result, irange, tmpl);
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return result;
}

SEXP _scan_bam(SEXP bfile, SEXP space, SEXP keepFlags, SEXP isSimpleCigar,
               SEXP reverseComplement, SEXP yieldSize, SEXP template_list,
               SEXP obeyQname, SEXP asMates)
{
    SEXP names = PROTECT(GET_ATTR(template_list, R_NamesSymbol));
    SEXP result =
        PROTECT(_scan_bam_result_init(template_list, names, space));
    SCAN_BAM_DATA sbd = _Calloc_SCAN_BAM_DATA(result);
    BAM_DATA bd = _init_BAM_DATA(bfile, space, keepFlags, isSimpleCigar,
                                 LOGICAL(reverseComplement)[0],
                                 INTEGER(yieldSize)[0],
                                 LOGICAL(obeyQname)[0], 
                                 LOGICAL(asMates)[0], (void *) sbd);

    int status = _do_scan_bam(bd, space, _filter_and_parse1,
                              _finish1range_BAM_DATA);
    if (status < 0) {
        int idx = bd->irec;
        int parse_status = bd->parse_status;
        _Free_SCAN_BAM_DATA(bd->extra);
        _Free_BAM_DATA(bd);
        Rf_error("'scanBam' failed:\n  record: %d\n  error: %d",
                 idx, parse_status);
    }

    _Free_SCAN_BAM_DATA(bd->extra);
    _Free_BAM_DATA(bd);
    UNPROTECT(2);
    return result;
}

/* count */

static int _count_bam1(const bam1_t * bam, void *data)
{
    return _count1_BAM_DATA(bam, (BAM_DATA) data);
}

SEXP _count_bam(SEXP bfile, SEXP space, SEXP keepFlags, SEXP isSimpleCigar)
{
    SEXP result = PROTECT(NEW_LIST(2));
    BAM_DATA bd =
        _init_BAM_DATA(bfile, space, keepFlags, isSimpleCigar, 0,
                       NA_INTEGER, 0, 0, result);

    SET_VECTOR_ELT(result, 0, NEW_INTEGER(bd->nrange));
    SET_VECTOR_ELT(result, 1, NEW_NUMERIC(bd->nrange));
    for (int i = 0; i < bd->nrange; ++i) {
        INTEGER(VECTOR_ELT(result, 0))[i] = REAL(VECTOR_ELT(result, 1))[i] = 0;
    }

    SEXP nms = PROTECT(NEW_CHARACTER(2));
    SET_STRING_ELT(nms, 0, mkChar("records"));
    SET_STRING_ELT(nms, 1, mkChar("nucleotides"));
    setAttrib(result, R_NamesSymbol, nms);
    UNPROTECT(1);

    int status = _do_scan_bam(bd, space, _count_bam1, NULL);
    if (status < 0) {
        int idx = bd->irec;
        int parse_status = bd->parse_status;
        _Free_BAM_DATA(bd);
        UNPROTECT(1);
        Rf_error("'countBam' failed:\n  record: %d\n  error: %d",
                 idx, parse_status);
    }

    _Free_BAM_DATA(bd);
    UNPROTECT(1);
    return result;
}

void scan_bam_cleanup()
{
    /* placeholder */
}

/* filterBam */

static int _prefilter_bam1(const bam1_t * bam, void *data)
{
    BAM_DATA bd = (BAM_DATA) data;
    bd->irec += 1;
    if (!_filter1_BAM_DATA(bam, bd))
        return 0;
    bambuffer_push((BAM_BUFFER) bd->extra, bam);
    bd->iparsed += 1;
    return 1;
}

SEXP
_prefilter_bam(SEXP bfile, SEXP space, SEXP keepFlags, SEXP isSimpleCigar,
               SEXP yieldSize, SEXP obeyQname, SEXP asMates)
{
    SEXP ext = PROTECT(bambuffer(INTEGER(yieldSize)[0]));
    BAM_DATA bd = _init_BAM_DATA(bfile, space, keepFlags, isSimpleCigar,
                                 0, INTEGER(yieldSize)[0],
                                 LOGICAL(obeyQname)[0], 
                                 LOGICAL(asMates), BAMBUFFER(ext));

    int status = _do_scan_bam(bd, space, _prefilter_bam1, NULL);
    if (status < 0) {
        int idx = bd->irec;
        int parse_status = bd->parse_status;
        _Free_BAM_DATA(bd);
        UNPROTECT(1);
        Rf_error("'filterBam' prefilter failed:\n  record: %d\n  error: %d",
                 idx, parse_status);
    }

    _Free_BAM_DATA(bd);
    UNPROTECT(1);
    return ext;
}

static int _filter_bam1(const bam1_t * bam, void *data)
{
    BAM_DATA bd = (BAM_DATA) data;
    bd->irec += 1;
    if (!_filter1_BAM_DATA(bam, bd))
        return 0;
    samwrite((samfile_t *) bd->extra, bam);
    bd->iparsed += 1;
    return 1;
}

SEXP
_filter_bam(SEXP bfile, SEXP space, SEXP keepFlags,
            SEXP isSimpleCigar, SEXP fout_name, SEXP fout_mode)
{
    /* open destination */
    BAM_DATA bd =
        _init_BAM_DATA(bfile, space, keepFlags, isSimpleCigar, 0,
                       NA_INTEGER, 0, 0, NULL);
    /* FIXME: this just copies the header... */
    bam_header_t *header = BAMFILE(bfile)->file->header;
    samfile_t *f_out = _bam_tryopen(translateChar(STRING_ELT(fout_name, 0)),
                                    CHAR(STRING_ELT(fout_mode, 0)), header);
    bd->extra = f_out;

    int status = _do_scan_bam(bd, space, _filter_bam1, NULL);
    if (status < 0) {
        int idx = bd->irec;
        int parse_status = bd->parse_status;
        _Free_BAM_DATA(bd);
        samclose(f_out);
        Rf_error("'filterBam' failed:\n  record: %d\n  error: %d",
                 idx, parse_status);
    }

    /* sort and index destintation ? */
    /* cleanup */
    _Free_BAM_DATA(bd);
    samclose(f_out);

    return status < 0 ? R_NilValue : fout_name;
}

/* merge_bam */

/* from bam_sort.c */

#define MERGE_RG     1
#define MERGE_LEVEL1 4
#define MERGE_FORCE  8

int bam_merge_core(int by_qname, const char *out, const char *headers,
                   int n, char * const *fn, int flag, const char *reg);

SEXP merge_bam(SEXP fnames, SEXP destination, SEXP overwrite,
               SEXP hname, SEXP regionStr, SEXP isByQname,
               SEXP addRG, SEXP compressLevel1)
{
    int i;

    if (!IS_CHARACTER(fnames) || 2 > Rf_length(fnames))
        Rf_error("'files' must be a character() with length >= 2");
    if (!IS_CHARACTER(hname) || 1 <  Rf_length(hname))
        Rf_error("'header' must be character() with length <= 1");
    if (!IS_CHARACTER(destination) || 1 != Rf_length(destination))
        Rf_error("'destination' must be character(1)");
    if (!IS_LOGICAL(overwrite) || 1 != Rf_length(overwrite))
        Rf_error("'overwrite' must be logical(1)");
    if (!IS_CHARACTER(regionStr) || 1 <  Rf_length(regionStr))
        Rf_error("'region' must define 0 or 1 regions");
    if (!IS_LOGICAL(isByQname) || 1 != Rf_length(isByQname))
        Rf_error("'isByQname' must be logical(1)");
    if (!IS_LOGICAL(addRG) || 1 != Rf_length(addRG))
        Rf_error("'addRG' must be logical(1)");
    if (!IS_LOGICAL(compressLevel1) || 1 != Rf_length(compressLevel1))
        Rf_error("'compressLevel1' must be logical(1)");

    char ** fileNames = (char **)
        R_alloc(sizeof(const char *), Rf_length(fnames));
    for (i = 0; i < Rf_length(fnames); ++i)
        fileNames[i] = (char *) translateChar(STRING_ELT(fnames, i));

    const char *hfName = 0 == Rf_length(hname) ?
        NULL : translateChar(STRING_ELT(hname, 0));

    int flag = 0;
    if (LOGICAL(addRG)[0])
        flag = ((int) flag) | ((int) MERGE_RG);
    if (LOGICAL(overwrite)[0])
        flag = ((int) flag) | ((int) MERGE_FORCE);
    if (LOGICAL(compressLevel1)[0])
        flag = ((int) flag) | ((int) MERGE_LEVEL1);

    const char *region = 0 == Rf_length(regionStr) ?
        NULL : translateChar(STRING_ELT(regionStr, 0));

    int res = bam_merge_core(LOGICAL(isByQname)[0],
                             translateChar(STRING_ELT(destination, 0)),
                             hfName, Rf_length(fnames), fileNames,
                             flag, region);
    if (res < 0)
        Rf_error("'mergeBam' failed with error code %d", res);

    return destination;
}

/* sort_bam */

SEXP sort_bam(SEXP filename, SEXP destination, SEXP isByQname, SEXP maxMemory)
{
    if (!IS_CHARACTER(filename) || 1 != LENGTH(filename))
        Rf_error("'filename' must be character(1)");
    if (!IS_CHARACTER(destination) || 1 != LENGTH(destination))
        Rf_error("'destination' must be character(1)");
    if (!IS_LOGICAL(isByQname) || LENGTH(isByQname) != 1)
        Rf_error("'isByQname' must be logical(1)");
    if (!IS_INTEGER(maxMemory) || LENGTH(maxMemory) != 1 ||
        INTEGER(maxMemory)[0] < 1)
        Rf_error("'maxMemory' must be a positive integer(1)");

    const char *fbam = translateChar(STRING_ELT(filename, 0));
    const char *fout = translateChar(STRING_ELT(destination, 0));
    int sortMode = asInteger(isByQname);

    size_t maxMem = (size_t) INTEGER(maxMemory)[0] * 1024 * 1024;
    _check_is_bam(fbam);
    bam_sort_core(sortMode, fbam, fout, maxMem);

    return destination;
}

/* index_bam */

SEXP index_bam(SEXP indexname)
{
    if (!IS_CHARACTER(indexname) || 1 != LENGTH(indexname))
        Rf_error("'indexname' must be character(1)");
    const char *fbam = translateChar(STRING_ELT(indexname, 0));

    _check_is_bam(fbam);
    int status = bam_index_build(fbam);

    if (0 != status)
        Rf_error("failed to build index\n  file: %s", fbam);
    char *fidx = (char *) R_alloc(strlen(fbam) + 5, sizeof(char));
    sprintf(fidx, "%s.bai", fbam);
    return mkString(fidx);
}
