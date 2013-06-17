// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file
///
/// \author Chris Saunders
///


#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "blt_util/parse_util.hh"
#include "blt_util/seq_util.hh"
#include "starling_common/align_path.hh"

#include "boost/foreach.hpp"
#include "boost/lexical_cast.hpp"

#include <cassert>
#include <cctype>

#include <algorithm>
#include <iostream>
#include <sstream>



enum {
    INDEL_BEGIN='^',
    INDEL_END='$'
};



static
void
unknown_md_error(const char* const md,
                 const char* const mdptr) {

    std::ostringstream oss;
    oss << "ERROR: can't parse match descriptor string: " << md << "\n"
        << "\tunexpected character: '" << *mdptr << "' at position: " << (mdptr-md+1) << "\n";
    throw blt_exception(oss.str().c_str());
}



static
void
unknown_cigar_error(const char* const cigar,
                    const char* const cptr) {

    std::ostringstream oss;
    oss << "ERROR: can't parse cigar string: " << cigar << "\n"
        << "\tunexpected character: '" << *cptr << "' at position: " << (cptr-cigar+1) << "\n";
    throw blt_exception(oss.str().c_str());
}



namespace ALIGNPATH {


static
void
apath_push(path_t& apath,
           path_segment& ps,
           const align_t t) {

    if( (0==ps.length) || (ps.type==t) ) return;
    apath.push_back(ps);
    ps.clear();
}



static
void
export_md_to_apath_impl(const char* md,
                        path_t& apath) {

    using illumina::blt_util::parse_unsigned;

    const char* mdptr(md);
    path_segment ps;

    while(*mdptr) {
        if       (isdigit(*mdptr)) {
            apath_push(apath,ps,MATCH);
            const unsigned mlen(parse_unsigned(mdptr));
            ps.length += mlen;
            ps.type = MATCH;

        } else if(is_valid_base(*mdptr)) {
            apath_push(apath,ps,MATCH);
            mdptr++;
            ps.length++;
            ps.type = MATCH;

        } else if(*mdptr == INDEL_BEGIN) {
            mdptr++; // eat INDEL_BEGIN

            while(*mdptr != INDEL_END) {
                if       (isdigit(*mdptr)) {
                    apath_push(apath,ps,INSERT);
                    const unsigned mlen(parse_unsigned(mdptr));
                    ps.length=mlen;
                    ps.type=INSERT;

                } else if(is_valid_base(*mdptr)) {
                    apath_push(apath,ps,DELETE);
                    mdptr++;
                    ps.length++;
                    ps.type=DELETE;

                } else {
                    unknown_md_error(md,mdptr);
                }
            }

            mdptr++; // eat INDEL_END

        } else {
            unknown_md_error(md,mdptr);
        }
    }

    apath_push(apath,ps,NONE);
}



void
export_md_to_apath(const char* md,
                   const bool is_fwd_strand,
                   path_t& apath,
                   const bool is_edge_deletion_error) {

    // to make best use of previous code, we parse the MD in the
    // alignment direction and then orient apath to the forward strand
    // as a second step if required
    //
    assert(NULL != md);

    apath.clear();
    export_md_to_apath_impl(md,apath);

    unsigned as(apath.size());

    if( ((as>0) and (apath.front().type == DELETE)) or
        ((as>1) and (apath.back().type == DELETE)) ) {
        std::ostringstream oss;
        if(is_edge_deletion_error) {
            oss << "ERROR: ";
        } else {
            oss << "WARNING: ";
        }
        oss << "alignment path: " << apath_to_cigar(apath) << " contains meaningless edge deletion.\n";
        if(is_edge_deletion_error) {
            throw blt_exception(oss.str().c_str());
        } else {
            log_os << oss.str();
            path_t apath2;
            for(unsigned i(0); i<as; ++i) {
                if(((i==0) or ((i+1)==as)) and
                   apath[i].type == DELETE) continue;
                apath2.push_back(apath[i]);
            }
            apath=apath2;
            as=apath.size();
        }
    }

    if( (not is_fwd_strand) and (as>1) ) {
        std::reverse(apath.begin(),apath.end());
    }
}



void
apath_to_cigar(const path_t& apath,
               std::string& cigar) {

    cigar.clear();

    const unsigned as(apath.size());
    for(unsigned i(0); i<as; ++i) {
        const path_segment& ps(apath[i]);
        cigar += boost::lexical_cast<std::string>(ps.length);
        cigar.push_back(segment_type_to_cigar_code(ps.type));
    }
}



static
void
fwd_apath_to_export_md(path_t& apath,
                       const char* ref_begin,
                       const char* ref_bases,
                       const char* ref_end,
                       const char* read_bases,
                       std::string& md) {

    // process the align path
    bool foundUnsupportedCigar = false;
    path_t::const_iterator pCIter;
    for(pCIter = apath.begin(); pCIter != apath.end(); ++pCIter) {

        if(pCIter->type == DELETE) {

            // handle deletion
            md.push_back('^');
            for(uint32_t i = 0; i < pCIter->length; ++i, ++ref_bases) {
                md.push_back(*ref_bases);
            }
            md.push_back('$');

        } else if(pCIter->type == INSERT) {

            // handle insertion
            md.push_back('^');
            md += boost::lexical_cast<std::string>(pCIter->length);
            read_bases += pCIter->length;
            md.push_back('$');

        } else if(pCIter->type == MATCH) {

            // handle match/mismatch
            uint32_t numMatchingBases = 0;
            for(uint32_t i = 0; i < pCIter->length; ++i, ++ref_bases, ++read_bases) {

                // handle circular genome
                if((ref_bases < ref_begin) || (ref_bases > ref_end)) {
                    md.push_back('N');
                    continue;
                }

                if(*ref_bases != *read_bases) {

                    // write the number of preceding matching bases
                    if(numMatchingBases != 0) {
                        md += boost::lexical_cast<std::string>(numMatchingBases);
                        numMatchingBases = 0;
                    }

                    // output the mismatched base
                    md.push_back(*ref_bases);

                } else ++numMatchingBases;
            }

            // write the number of trailing matching bases
            if(numMatchingBases != 0) {
                md += boost::lexical_cast<std::string>(numMatchingBases);
            }

        } else {

            // handle unsupported CIGAR operation
            foundUnsupportedCigar = true;
            break;
        }
    }

    if(foundUnsupportedCigar) md = "UNSUPPORTED";
}



static
void
rev_apath_to_export_md(path_t& apath,
                       const char* ref_begin,
                       const char* ref_bases,
                       const char* ref_end,
                       const char* read_bases,
                       std::string& md) {

    // process the align path
    bool foundUnsupportedCigar = false;
    path_t::const_reverse_iterator pCRIter;
    for(pCRIter = apath.rbegin(); pCRIter != apath.rend(); ++pCRIter) {

        if(pCRIter->type == DELETE) {

            // handle deletion
            md.push_back('^');
            for(uint32_t i = 0; i < pCRIter->length; ++i, --ref_bases) {
                md.push_back(comp_base(*ref_bases));
            }
            md.push_back('$');

        } else if(pCRIter->type == INSERT) {

            // handle insertion
            md.push_back('^');
            md += boost::lexical_cast<std::string>(pCRIter->length);
            read_bases += pCRIter->length;
            md.push_back('$');

        } else if(pCRIter->type == MATCH) {

            // recreate the the match descriptor for this non-INDEL region
            uint32_t numMatchingBases = 0;
            for(uint32_t i = 0; i < pCRIter->length; ++i, --ref_bases, ++read_bases) {

                // handle circular genome
                if((ref_bases < ref_begin) || (ref_bases > ref_end)) {
                    md.push_back('N');
                    continue;
                }

                const char rcRefBase = comp_base(*ref_bases);

                if(rcRefBase != *read_bases) {

                    // write the number of preceding matching bases
                    if(numMatchingBases != 0) {
                        md += boost::lexical_cast<std::string>(numMatchingBases);
                        numMatchingBases = 0;
                    }

                    // output the mismatched base
                    md.push_back(rcRefBase);

                } else ++numMatchingBases;
            }

            // write the number of trailing matching bases
            if(numMatchingBases != 0) {
                md += boost::lexical_cast<std::string>(numMatchingBases);
            }

        } else {

            // handle unsupported CIGAR operation
            foundUnsupportedCigar = true;
            break;
        }
    }

    if(foundUnsupportedCigar) md = "UNSUPPORTED";
}



void
apath_to_export_md(path_t& apath,
                   const char* ref_seq,
                   const char* ref_end,
                   const int32_t ref_pos,
                   const std::string& read_bases,
                   const bool is_fwd_strand,
                   std::string& md) {

    md.clear();

    if(is_fwd_strand) {

        const char* pRead      = read_bases.c_str();
        const char* pReference = ref_seq + ref_pos - 1;
        fwd_apath_to_export_md(apath, ref_seq, pReference, ref_end, pRead, md);

    } else {

        uint32_t numRefBases = 0;
        path_t::const_iterator pCIter;
        for(pCIter = apath.begin(); pCIter != apath.end(); ++pCIter) {
            if((pCIter->type == MATCH) || (pCIter->type == DELETE)) {
                numRefBases += pCIter->length;
            }
        }

        const char* pRead      = read_bases.c_str();
        const char* pReference = ref_seq + ref_pos + numRefBases - 2;
        rev_apath_to_export_md(apath, ref_seq, pReference, ref_end, pRead, md);
    }
}



std::ostream&
operator<<(std::ostream& os, const path_t& apath) {
    const unsigned as(apath.size());
    for(unsigned i(0); i<as; ++i) {
        os << apath[i].length << segment_type_to_cigar_code(apath[i].type);
    }
    return os;
}



void
cigar_to_apath(const char* cigar,
               path_t& apath) {

    using illumina::blt_util::parse_unsigned;

    assert(NULL != cigar);

    apath.clear();

    path_segment lps;
    const char* cptr(cigar);
    while(*cptr) {
        path_segment ps;
        // expect sequences of digits and cigar codes:
        if(! isdigit(*cptr)) unknown_cigar_error(cigar,cptr);
        ps.length = parse_unsigned(cptr);
        ps.type = cigar_code_to_segment_type(*cptr);
        if(ps.type == NONE) unknown_cigar_error(cigar,cptr);
        cptr++;
        if((ps.type == PAD) || (ps.length == 0)) continue;

        if(ps.type != lps.type) {
            if(lps.type != NONE) apath.push_back(lps);
            lps = ps;
        } else {
            lps.length += ps.length;
        }
    }

    if(lps.type != NONE) apath.push_back(lps);
}



unsigned
apath_read_length(const path_t& apath) {

    unsigned val(0);
    const unsigned as(apath.size());
    for(unsigned i(0); i<as; ++i) {
        const path_segment& ps(apath[i]);
        if(! is_segment_type_read_length(ps.type)) continue;
        val += ps.length;
    }
    return val;
}



unsigned
apath_ref_length(const path_t& apath) {

    unsigned val(0);
    const unsigned as(apath.size());
    for(unsigned i(0); i<as; ++i) {
        const path_segment& ps(apath[i]);
        if(! is_segment_type_ref_length(ps.type)) continue;
        val += ps.length;
    }
    return val;
}



static
inline
bool
is_segment_type_unaligned_read_edge(const align_t id) {
    switch(id) {
    case INSERT    :
    case HARD_CLIP :
    case SOFT_CLIP : return true;
    default        : return false;
    }
}



unsigned
apath_read_lead_size(const path_t& apath) {
    unsigned val(0);
    const unsigned as(apath.size());
    for(unsigned i(0); i<as; ++i) {
        const path_segment& ps(apath[i]);
        if(! is_segment_type_unaligned_read_edge(ps.type)) return val;
        if(is_segment_type_read_length(ps.type)) val += ps.length;
    }
    return val;
}



unsigned
apath_read_trail_size(const path_t& apath) {
    unsigned val(0);
    const unsigned as(apath.size());
    for(unsigned i(0); i<as; ++i) {
        const path_segment& ps(apath[as-i-1]);
        if(! is_segment_type_unaligned_read_edge(ps.type)) return val;
        if(is_segment_type_read_length(ps.type)) val += ps.length;
    }
    return val;
}



unsigned
apath_soft_clip_lead_size(const path_t& apath) {
    unsigned val(0);
    const unsigned as(apath.size());
    for(unsigned i(0); i<as; ++i) {
        const path_segment& ps(apath[i]);
        if       (HARD_CLIP == ps.type) {
            // do nothing:
        } else if(SOFT_CLIP == ps.type) {
            val += ps.length;
        } else {
            break;
        }
    }
    return val;
}



unsigned
apath_soft_clip_trail_size(const path_t& apath) {
    unsigned val(0);
    const unsigned as(apath.size());
    for(unsigned i(0); i<as; ++i) {
        const path_segment& ps(apath[as-i-1]);
        if       (HARD_CLIP == ps.type) {
            // do nothing:
        } else if(SOFT_CLIP == ps.type) {
            val += ps.length;
        } else {
            break;
        }
    }
    return val;
}



unsigned
apath_insert_lead_size(const path_t& apath) {
    unsigned val(0);
    const unsigned as(apath.size());
    for(unsigned i(0); i<as; ++i) {
        const path_segment& ps(apath[i]);
        if((HARD_CLIP == ps.type) || (SOFT_CLIP == ps.type)) {
            // do nothing:
        } else if(INSERT == ps.type) {
            val += ps.length;
        } else {
            break;
        }
    }
    return val;
}



unsigned
apath_insert_trail_size(const path_t& apath) {
    unsigned val(0);
    const unsigned as(apath.size());
    for(unsigned i(0); i<as; ++i) {
        const path_segment& ps(apath[as-i-1]);
        if((HARD_CLIP == ps.type) || (SOFT_CLIP == ps.type)) {
            // do nothing:
        } else if(INSERT == ps.type) {
            val += ps.length;
        } else {
            break;
        }
    }
    return val;
}



void
apath_clip_clipper(path_t& apath,
                   unsigned& hc_lead,
                   unsigned& hc_trail,
                   unsigned& sc_lead,
                   unsigned& sc_trail) {

    hc_lead=0;
    hc_trail=0;
    sc_lead=0;
    sc_trail=0;

    bool is_lead(true);
    path_t apath2;
    const unsigned as(apath.size());
    for(unsigned i(0); i<as; ++i) {
        const path_segment& ps(apath[i]);
        if       (HARD_CLIP == ps.type) {
            if(is_lead) { hc_lead += ps.length; }
            else        { hc_trail += ps.length; }
        } else if(SOFT_CLIP == ps.type) {
            if(is_lead) { sc_lead += ps.length; }
            else        { sc_trail += ps.length; }
        } else {
            is_lead=false;
            assert(0==hc_trail);
            assert(0==sc_trail);
            apath2.push_back(ps);
        }
    }
    apath=apath2;
}



void
apath_clip_adder(path_t& apath,
                 const unsigned hc_lead,
                 const unsigned hc_trail,
                 const unsigned sc_lead,
                 const unsigned sc_trail) {

    path_t apath2;
    path_segment ps;
    if(hc_lead>0) {
        ps.type = HARD_CLIP;
        ps.length = hc_lead;
        apath2.push_back(ps);
    }
    if(sc_lead>0) {
        ps.type = SOFT_CLIP;
        ps.length = sc_lead;
        apath2.push_back(ps);
    }
    apath2.insert(apath2.end(),apath.begin(),apath.end());
    if(sc_trail>0) {
        ps.type = SOFT_CLIP;
        ps.length = sc_trail;
        apath2.push_back(ps);
    }
    if(hc_trail>0) {
        ps.type = HARD_CLIP;
        ps.length = hc_trail;
        apath2.push_back(ps);
    }
    apath=apath2;
}



// 1. remove zero-length segments
// 2. remove pads
// 3. condense repeated segment types
// 4. reduce adjacent insertion/deletion tags to a single pair
//
// return true if path has been altered
//
bool
apath_cleaner(path_t& apath) {
    bool is_cleaned(false);
    const unsigned as(apath.size());
    unsigned insertIndex(as);
    unsigned deleteIndex(as);
    unsigned otherIndex(as);
    for(unsigned i(0); i<as; ++i) {
        path_segment& ps(apath[i]);
        if       (ps.length == 0) {
            is_cleaned = true;
        } else if(ps.type == PAD) {
            ps.length = 0;
            is_cleaned = true;
        } else if(ps.type == INSERT) {
            if(insertIndex < as) {
                apath[insertIndex].length += ps.length;
                ps.length = 0;
                is_cleaned = true;
            } else {
                insertIndex = i;
            }
        } else if(ps.type == DELETE) {
            if(deleteIndex < as) {
                apath[deleteIndex].length += ps.length;
                ps.length = 0;
                is_cleaned = true;
            } else {
                deleteIndex = i;
            }
        } else {
            if((insertIndex<as) || (deleteIndex<as)) {
                insertIndex = as;
                deleteIndex = as;
                otherIndex = as;
            }
            if((otherIndex < as) && (apath[otherIndex].type == ps.type)) {
                apath[otherIndex].length += ps.length;
                ps.length = 0;
                is_cleaned = true;
            } else {
                otherIndex = i;
            }
        }
    }
    if(is_cleaned) {
        path_t apath2;
        for(unsigned i(0); i<as; ++i) {
            if(apath[i].length == 0) continue;
            apath2.push_back(apath[i]);
        }
        apath = apath2;
    }
    return is_cleaned;
}


#if 0
std::pair<unsigned,unsigned>
get_nonclip_end_segments(const path_t& apath) {
    const unsigned as(apath.size());
    std::pair<unsigned,unsigned> res(as,as);
    bool is_first_nonclip(false);
    for(unsigned i(0); i<as; ++i) {
        const path_segment& ps(apath[i]);
        if(! (ps.type == SOFT_CLIP ||
              ps.type == HARD_CLIP)) {
            if(! is_first_nonclip) {
                res.first=i;
                is_first_nonclip=true;
            }
            res.second=i;
        }
    }
    return res;
}
#endif


pos_range
get_nonclip_range(const path_t& apath) {
    pos_range pr;
    unsigned read_offset(0);
    const unsigned as(apath.size());
    for(unsigned i(0); i<as; ++i) {
        const path_segment& ps(apath[i]);
        const bool is_rt(is_segment_type_read_length(ps.type));
        if(! (ps.type == SOFT_CLIP ||
              ps.type == HARD_CLIP)) {
            if(! pr.is_begin_pos) {
                pr.set_begin_pos(read_offset);
            }
            pr.set_end_pos(read_offset + (is_rt ? ps.length : 0));
        }
        if(is_rt) read_offset+=ps.length;
    }
    return pr;
}



std::pair<unsigned,unsigned>
get_match_edge_segments(const path_t& apath) {
    const unsigned as(apath.size());
    std::pair<unsigned,unsigned> res(as,as);
    bool is_first_match(false);
    for(unsigned i(0); i<as; ++i) {
        const path_segment& ps(apath[i]);
        if(MATCH == ps.type) {
            if(! is_first_match) res.first=i;
            is_first_match=true;
            res.second=i;
        }
    }
    return res;
}



unsigned
apath_exon_count(const path_t& apath) {
    unsigned val(1);
    const unsigned as(apath.size());
    for(unsigned i(0); i<as; ++i) {
        if(apath[i].type==SKIP) val++;
    }
    return val;
}



bool
is_clipped(const path_t& apath) {
    const unsigned as(apath.size());
    if(as==0) return false;
    if((apath[0].type == SOFT_CLIP) || (apath[0].type == HARD_CLIP)) return true;
    if(as>1) {
        if((apath[as-1].type == SOFT_CLIP) || (apath[as-1].type == HARD_CLIP)) return true;
    }
    return false;
}



bool
is_soft_clipped(const path_t& apath) {
    BOOST_FOREACH(const path_segment& ps, apath) {
        if(SOFT_CLIP == ps.type) return true;
    }
    return false;
}



bool
is_edge_readref_len_segment(const path_t& apath) {
    const unsigned as(apath.size());
    if(as==0) return false;

    const std::pair<unsigned,unsigned> ends(get_match_edge_segments(apath));

    // at this point we assume the alignment has been sanity checked for legal clipping,
    // where hard-clip is only on the outside, next soft-clipping, then anything else...
    //
    for(unsigned i(0); i<as; ++i) {
        const path_segment& ps(apath[i]);

        const bool is_edge_segment((i<ends.first) || (i>ends.second));
        const bool is_clip_type(ps.type==INSERT || ps.type==DELETE || ps.type==SKIP || ps.type==SOFT_CLIP);
        if(is_edge_segment && is_clip_type) return true;
    }
    return false;
}



bool
is_seq_swap(const path_t& apath) {
    const unsigned as(apath.size());
    for(unsigned i(0); (i+1)<as; ++i) {
        if(is_segment_type_indel(apath[i].type) &&
           is_segment_type_indel(apath[i+1].type)) {
            return true;
        }
    }
    return false;
}



bool
is_segment_swap_start(const path_t& apath,
                      unsigned i) {

    using namespace ALIGNPATH;

    bool is_insert(false);
    bool is_delete(false);

    const unsigned as(apath.size());
    for(;i<as;++i) {
        if     (apath[i].type == INSERT) { is_insert=true; }
        else if(apath[i].type == DELETE) { is_delete=true; }
        else { break; }
    }

    return (is_insert && is_delete);
}



bool
is_apath_floating(const path_t& apath) {

    BOOST_FOREACH(const path_segment& ps, apath) {
        if(ps.type==MATCH) return false;
    }
    return true;
}


std::string
get_apath_invalid_reason(const path_t& apath,
                         const unsigned seq_length) {

    const ALIGN_ISSUE::issue_t ai(get_apath_invalid_type(apath,seq_length));

    if(ALIGN_ISSUE::LENGTH == ai) {
        std::ostringstream oss;
        oss << "alignment length (" << apath_read_length(apath) << ") does not match read length (" << seq_length << ")";
        return oss.str();
    }

    return std::string(ALIGN_ISSUE::description(ai));
}



ALIGN_ISSUE::issue_t
get_apath_invalid_type(const path_t& apath,
                     const unsigned seq_length) {

    bool is_match(false);
    align_t last_type(NONE);
    const unsigned as(apath.size());
    for(unsigned i(0); i<as; ++i) {
        const path_segment& ps(apath[i]);

        if(ps.type==NONE) return ALIGN_ISSUE::UNKNOWN_SEGMENT;
        if((i!=0) && ps.type==last_type) return ALIGN_ISSUE::REPEATED_SEGMENT;

        if(! is_match) {
            if(ps.type==SKIP) return ALIGN_ISSUE::EDGE_SKIP;
        }

        if(ps.type==HARD_CLIP) {
            if(! ((i==0) || ((i+1)==as))) return ALIGN_ISSUE::CLIPPING;
        }

        if(ps.type==SOFT_CLIP) {
            if(! ((i==0) || ((i+1)==as))) {
                if(i==1) {
                    if(as==3) {
                        if((apath[0].type != HARD_CLIP) && (apath[i+1].type != HARD_CLIP)) return ALIGN_ISSUE::CLIPPING;
                    } else {
                        if(apath[0].type != HARD_CLIP) return ALIGN_ISSUE::CLIPPING;
                    }
                } else if((i+2)==as) {
                    if(apath[i+1].type != HARD_CLIP) return ALIGN_ISSUE::CLIPPING;
                } else {
                    return ALIGN_ISSUE::CLIPPING;
                }
            }
        }

        if((! is_match) && (ps.type==MATCH)) is_match=true;

        last_type=ps.type;
    }

    if(! is_match) return ALIGN_ISSUE::FLOATING;

    // run in reverse to finish checking condition (2a):
    for(unsigned i(0); i<as; ++i) {
        const path_segment& ps(apath[as-(i+1)]);
        if(ps.type==MATCH) break;
        //if(ps.type==DELETE) return ALIGN_ISSUE::EDGE_DELETE;
        if(ps.type==SKIP) return ALIGN_ISSUE::EDGE_SKIP;
    }

    if(seq_length != apath_read_length(apath)) return ALIGN_ISSUE::LENGTH;

    return ALIGN_ISSUE::NONE;
}



// Unlike the above function which tests for invalid alignment paths,
// this function test for valid alignment methods which starling
// simply cannot handle
//
bool
is_apath_starling_invalid(const path_t& apath) {

    BOOST_FOREACH(const path_segment& ps, apath) {
        if(ps.type==PAD) return true;
    }
    return false;
}



}  // namespace ALIGNPATH
