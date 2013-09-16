// BamFileIterator.hpp:
// Iterator used when reading ranges from a bam file.

#ifndef BAMRANGEITERATOR_H
#define BAMRANGEITERATOR_H

#include "BamIterator.hpp"


class BamRangeIterator : public BamIterator {

    bam_iter_t iter;

    void iterate_complete(bamFile bfile) {
	if (iter_done)
	    return;
	if (NULL == bam) {	// first record 
	    bam = bam_init1();
	    if (bam_iter_read(bfile, iter, bam) < 0) {
		iter_done = true;
		return;
	    }
	}

	bool done = false;
	do {
	    process(bam);
            int tid = bam->core.tid;
            int pos = bam->core.pos;
	    if (bam_iter_read(bfile, iter, bam) < 0) {
		iter_done = done = true;
            } else if (complete.size() != 0) {
		// stop if something to yield AND finished position
		    done = (bam->core.tid != tid) || (bam->core.pos != pos);
		}
	} while (!done);
    }

    void iterate_incomplete(bamFile bfile) {
        int64_t pos = bam_tell(bfile);
        Templates::iterator it;
        for (it = templates.begin(); it != templates.end(); ++it) {
            // mate all segments in 'inprogress'
            while (it->second.mate_inprogress_segments(bfile, bindex))
                complete.push(it->second.get_complete());
        }
        BamIterator::iterate_incomplete(bfile);
        bam_seek(bfile, pos, SEEK_SET);
    }

public:

    // constructor / destructor
    BamRangeIterator(const bam_index_t *bindex,
                     int tid, int beg, int end) :
        BamIterator(bindex) 
    {
	iter = bam_iter_query(bindex, tid, beg, end);
    }

    ~BamRangeIterator() {
	bam_iter_destroy(iter);
   }
};

#endif
