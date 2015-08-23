// BamRangeIterator.h:
// Iterator used when reading ranges from a bam file.

#ifndef BAMRANGEITERATOR_H
#define BAMRANGEITERATOR_H

#include "BamIterator.h"

class BamRangeIterator : public BamIterator {

    int32_t tid, beg, end;
    hts_itr_t *iter;

    void iterate_inprogress(samFile *bfile) {
	if (NULL == bam) {	// first record 
	    bam = bam_init1();
	    if (sam_itr_next(bfile, iter, bam) < 0) {
		iter_done = true;
		return;
	    }
	}

	do {
	    process(bam);
	    if (sam_itr_next(bfile, iter, bam) < 0)
		iter_done = true;
	} while (!iter_done);
        mate_touched_templates();
    }

    void finalize_inprogress(samFile *bfile) {
        int64_t pos = bgzf_tell(bfile->fp.bgzf);
        Templates::iterator it;
        // mate 'inprogress' segments for all templates
        for (it = templates.begin(); it != templates.end(); ++it)
            it->second.mate_inprogress_segments(bfile, bindex, complete,
                                                qname_prefix_end(),
                                                qname_suffix_start(),
                                                tid, beg, end,
                                                header->target_len,
                                                it->first);

        BamIterator::finalize_inprogress(bfile);
        bgzf_seek(bfile->fp.bgzf, pos, SEEK_SET);
    }

public:

    // constructor / destructor
    BamRangeIterator(samFile *bfile, const hts_idx_t *bindex,
                     int32_t tid, int32_t beg, int32_t end) :
        BamIterator(bfile, bindex), tid(tid), beg(beg), end(end)
    {
	iter = bam_itr_queryi(bindex, tid, beg == 0 ? 0 : beg - 1, end);
    }

    ~BamRangeIterator() {
	sam_itr_destroy(iter);
   }
};

#endif
