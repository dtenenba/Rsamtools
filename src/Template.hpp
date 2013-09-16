#ifndef TEMPLATE_H
#define TEMPALTE_H

#include <algorithm>
#include <list>
#include "samtools/sam.h"

using namespace std;

class Template {
    typedef list<const bam1_t *> Segments;
    typedef Segments::iterator iterator;
    typedef Segments::const_iterator const_iterator;

    char *rg, *qname;
    Segments inprogress, complete, incomplete; 

    // FIXME: check RG retrieval
    const int readgroup_q(const bam1_t *mate) const {
	uint8_t *aux = bam_aux_get(mate, "RG");
	char *rg0 = '\0';
	if (aux != 0)
	    rg0 = bam_aux2Z(aux);

        if (rg =='\0' && rg0 == '\0') /* both null ok */
            return 0;
        else
	return strcmp(rg, rg0);
    }

    const int qname_q(const bam1_t *mate) const {
	return strcmp(qname, bam1_qname(mate));
    }

public:

    Template() : rg('\0') {}

    Template(const bam1_t *bam) : 
        rg('\0') { add_segment(bam); }

    // getters / setters
    size_t size() const { 
        return (inprogress.size() + complete.size() + incomplete.size()); 
    }

    list<const bam1_t *> get_incomplete() {
        Segments segments(incomplete);
        incomplete.clear();
        return segments;
    }

    list<const bam1_t *> get_complete() {
        Segments segments(complete);
        complete.clear();
        return segments;
    }

    // is_valid
    // 1. Bit 0x1 (multiple segments) is 1 
    // 2. Bit 0x4 (segment unmapped) is 0
    // 3. Bit 0x8 (next segment unmapped) is 0
    // 4. 'mpos' != -1 (i.e. PNEXT = 0)
    bool is_valid(const bam1_t *bam) {
        const bool multi_seg = bam->core.flag & BAM_FPAIRED;
        const bool seg_unmapped = bam->core.flag & BAM_FUNMAP;
        const bool mate_unmapped = bam->core.flag & BAM_FMUNMAP;

        return  multi_seg && !seg_unmapped && 
                !mate_unmapped && (bam->core.mpos != -1);
    }

    // is_template
    bool is_template(const bam1_t *mate) const {
	return (readgroup_q(mate) == 0) && (qname_q(mate) == 0);
    }

    // is_mate
    // 1. Bit 0x40 and 0x80: Segments are a pair of first/last OR
    //    neither segment is marked first/last
    // 2. Bit 0x100: Both segments are secondary OR both not secondary
    // 3. tid match
    // 4. segment1 mpos matches segment2 pos AND
    //    segment2 mpos matches segment1 pos
    bool is_mate(const bam1_t *bam, const bam1_t *mate) const {
        const bool bam_read1 = bam->core.flag & BAM_FREAD1;
        const bool bam_read2 = bam->core.flag & BAM_FREAD2;
        const bool bam_secondary = bam->core.flag & BAM_FSECONDARY;
        const bool mate_read1 = mate->core.flag & BAM_FREAD1;
        const bool mate_read2 = mate->core.flag & BAM_FREAD2;
        const bool mate_secondary = mate->core.flag & BAM_FSECONDARY;
	return
            ((bam_read1 == mate_read2) || (bam_read2 == mate_read1)) &&
            (bam_secondary == mate_secondary) &&
	    (bam->core.mtid == mate->core.tid) &&
            (bam->core.pos == mate->core.mpos) &&
	    (bam->core.mpos == mate->core.pos);
    }

    // add_segment
    // Returns true if segment added was a mate
    bool add_segment(const bam1_t *bam1) {
        bam1_t *bam = bam_dup1(bam1);
        if (!is_valid(bam)) {
            incomplete.push_back(bam);
	    return false;
        }

        // new template
        if (size() == 0) {
            qname = bam1_qname(bam);
            uint8_t *aux = bam_aux_get(bam, "RG");
            if (aux != 0)
                rg = bam_aux2Z(aux);
            inprogress.push_back(bam);
        // existing template with 'inprogress' records
        } else if (inprogress.size() > 0) {
            Segments::iterator it = inprogress.begin();
            while (it != inprogress.end()) {
                if (is_mate(bam, *it)) {
                    complete.push_back(*it);
                    complete.push_back(bam);
                    inprogress.erase(it); 
                    return true;
                } else it++;
            }
            inprogress.push_back(bam);
        // existing template with no 'inprogress' records
        } else {
            inprogress.push_back(bam);
        }
        return false;
    }

    // mate (BamRangeIterator only)
    bool mate_inprogress_segments(bamFile bfile, const bam_index_t * bindex) {
        // search for mate for each 'inprogress' segment
        bam1_t *bam = bam_init1();

        iterator it = inprogress.begin();
        while (it != inprogress.end()) {
            bool mate_found = false;
            const bam1_t *curr = *it;
            const int tid = curr->core.mtid;
            const int beg = curr->core.mpos;

            bam_iter_t iter = bam_iter_query(bindex, tid, beg, beg + 1);
            while (bam_iter_read(bfile, iter, bam) >= 0) {
                if (beg == -1)
                    break;      // no mate available
                if (is_valid(bam) && is_template(bam) && is_mate(curr, bam)) {
                    mate_found = true;
                    break;
                }
            }
            bam_iter_destroy(iter);
            if (mate_found) {
                iterator it0 = it++;
                complete.push_back(bam_dup1(bam));
                complete.push_back(curr);
                inprogress.erase(it0);
            } else it++;
        }

        bam_destroy1(bam);
        return complete.size() > 0;
    }

    // cleanup
    // move 'inprogress' to 'incomplete', return 'incomplete'
    list<const bam1_t *> cleanup() {
        if (complete.size())
            Rf_error("Error in cleanup: 'complete' not empty");

        if  (!inprogress.empty()) {
            incomplete.splice(incomplete.end(), inprogress);
            inprogress.clear();
        }

        return get_incomplete();
    }

};

#endif
