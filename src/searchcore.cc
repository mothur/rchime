/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2026, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
  All rights reserved.

  Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
  Department of Informatics, University of Oslo,
  PO Box 1080 Blindern, NO-0316 Oslo, Norway

  This software is dual-licensed and available under a choice
  of one of two licenses, either under the terms of the GNU
  General Public License version 3 or the BSD 2-Clause License.


  GNU General Public License version 3

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.


  The BSD 2-Clause License

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  1. Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
  COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.

*/

#include "chime.h"
#include "vsearch.h"
#include "align_simd.h"
#include "dbindex.h"
#include "linmemalign.h"
#include "minheap.h"
#include "unique.h"
#include "seqcmp.h"
#include "searchcore.h"
#include "cpu.h"
#include <algorithm>  // std::count_if, std::min, std::max
#include <array>
#include <cassert>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cmath>  // std::pow
#include <cstdint> // int64_t, uint64_t
#include <cstdio>  // std::sscanf, std::size_t
#include <cstdlib>  // std::qsort
#include <cstring>  // std::strlen, std::memset, std::strcmp
#include <limits>
#include <vector>


/* per thread data */
/******************************************************************************/
inline auto hit_compare_byid_typed(struct hit const * lhs, struct hit const * rhs) -> int
{
    /*
     Order:
     accepted, then rejected (weak)
     high id, then low id
     early target, then late target
     */

    if (lhs->rejected < rhs->rejected)
    {
        return -1;
    }
    if (lhs->rejected > rhs->rejected)
    {
        return +1;
    }
    if (lhs->aligned > rhs->aligned)
    {
        return -1;
    }
    if (lhs->aligned < rhs->aligned)
    {
        return +1;
    }
    if (lhs->aligned == 0)
    {
        return 0;
    }
    if (lhs->id > rhs->id)
    {
        return -1;
    }
    if (lhs->id < rhs->id)
    {
        return +1;
    }
    if (lhs->target < rhs->target)
    {
        return -1;
    }
    if (lhs->target > rhs->target)
    {
        return +1;
    }
    return 0;
}
/******************************************************************************/
auto hit_compare_byid(const void * lhs, const void * rhs) -> int
{
    return hit_compare_byid_typed((struct hit *) lhs, (struct hit *) rhs);
}

/******************************************************************************/
Vsearch_Searchcore::Vsearch_Searchcore() {
    opts = Vsearch_Options::getInstance();
}
/******************************************************************************/
Vsearch_Searchcore::~Vsearch_Searchcore() {}
/******************************************************************************/
auto Vsearch_Searchcore::make_hits_span(struct searchinfo_s const * search_info) -> Span<struct hit> {
    assert(search_info != nullptr);
    assert(search_info->hit_count >= 0);
    auto const length = static_cast<std::size_t>(search_info->hit_count);
    return Span<struct hit>{search_info->hits, length};
}
/******************************************************************************/
auto Vsearch_Searchcore::count_number_of_hits_to_keep(struct searchinfo_s const * search_info) -> std::size_t {
    if (search_info == nullptr) { return std::size_t{0}; }
    auto const hits = make_hits_span(search_info);
    return static_cast<std::size_t>(std::count_if(hits.cbegin(), hits.cend(),
                                                  [](struct hit const & hit) -> bool {
                                                    return hit.accepted or hit.weak;
                                                  }));
}
/******************************************************************************/
auto Vsearch_Searchcore::copy_over_hits_to_be_kept(std::vector<struct hit> & hits,
                                 struct searchinfo_s const * search_info) -> void {
    if (search_info == nullptr) { return; }
    for (auto const & hit : make_hits_span(search_info)) {
      if (hit.accepted or hit.weak) {
        hits.emplace_back(hit);
      }
    }
}
/******************************************************************************/
auto Vsearch_Searchcore::free_rejected_alignments(struct searchinfo_s const * search_info) -> void {
    if (search_info == nullptr) { return; }
    for (auto const & hit : make_hits_span(search_info)) {
      if (not (hit.accepted or hit.weak) and hit.aligned) {
        util.xfree(hit.nwalignment);
      }
    }
}
/******************************************************************************/
auto Vsearch_Searchcore::topscores(struct searchinfo_s * searchinfo,
                                   Vsearch_Database* db,
                                   Vsearch_DBIndex* dbindex,
                                   Vsearch_Minheap* minheap) -> void
{
  /*
    Count the kmer hits in each database sequence and
    make a sorted list of a given number (th)
    of the database sequences with the highest number of matching kmers.
    These are stored in the min heap array.
  */

  /* count kmer hits in the database sequences */
  const int indexed_count = dbindex->getcount();

  /* zero counts */
  std::memset(searchinfo->kmers, 0, indexed_count * sizeof(count_t));

  minheap->clear(searchinfo->m);

  for (auto i = 0U; i < searchinfo->kmersamplecount; i++)
    {
      auto const kmer = searchinfo->kmersample[i];
      auto * bitmap = dbindex->getbitmap(kmer);

  if (bitmap != nullptr)
        {
#ifdef __x86_64__
          if (ssse3_present != 0)
            {
              cpu.increment_counters_from_bitmap_ssse3(searchinfo->kmers,
                                                   bitmap, indexed_count);
            }
          else
            {
              cpu.increment_counters_from_bitmap_sse2(searchinfo->kmers,
                                                  bitmap, indexed_count);
            }
#else
          cpu.increment_counters_from_bitmap(searchinfo->kmers, bitmap, indexed_count);
#endif
        }
      else
        {
          auto * list = dbindex->getmatchlist(kmer);
          auto const count = dbindex->getmatchcount(kmer);
          for (auto j = 0U; j < count; j++)
            {
              searchinfo->kmers[list[j]]++;
            }
        }
    }

  auto const minmatches = std::min(static_cast<unsigned int>(opts->opt_minwordmatches), searchinfo->kmersamplecount);

    for (auto i = 0; i < indexed_count; i++)
    {
      auto const count = searchinfo->kmers[i];

      if (count >= minmatches)
        {
          auto const seqno = dbindex->getmapping(i);
          unsigned int const length = db->getsequencelen(seqno);

          elem_t novel;
          novel.count = count;
          novel.seqno = seqno;
          novel.length = length;

          minheap->add(searchinfo->m, & novel);
        }
    }
  minheap->sort(searchinfo->m);
}

/******************************************************************************/
auto Vsearch_Searchcore::align_trim(struct hit * hit) -> void
{
  /* trim alignment and fill in info */
  /* assumes that the hit has been aligned */

  /* info for semi-global alignment (without gaps at ends) */

  hit->trim_aln_left = 0;
  hit->trim_q_left = 0;
  hit->trim_t_left = 0;
  hit->trim_aln_right = 0;
  hit->trim_q_right = 0;
  hit->trim_t_right = 0;

  /* left trim alignment */

  auto * p = hit->nwalignment;
  auto op = '\0';
  int64_t run = 0;
  if (*p != 0)
    {
      run = 1;
      auto scanlength = 0;
      sscanf(p, "%" PRId64 "%n", &run, &scanlength);
      op = *(p + scanlength);
      if (op != 'M')
        {
          hit->trim_aln_left = 1 + scanlength;
          if (op == 'D')
            {
              hit->trim_q_left = run;
            }
          else
            {
              hit->trim_t_left = run;
            }
        }
    }

  /* right trim alignment */

  auto * e = hit->nwalignment + strlen(hit->nwalignment);
  if (e > hit->nwalignment)
    {
      p = e - 1;
      op = *p;
      if (op != 'M')
        {
          while ((p > hit->nwalignment) and (*(p - 1) <= '9'))
            {
              --p;
            }
          run = 1;
          sscanf(p, "%" PRId64, &run);
          hit->trim_aln_right = e - p;
          if (op == 'D')
            {
              hit->trim_q_right = run;
            }
          else
            {
              hit->trim_t_right = run;
            }
        }
    }

  if (hit->trim_q_left >= hit->nwalignmentlength)
    {
      hit->trim_q_right = 0;
    }

  if (hit->trim_t_left >= hit->nwalignmentlength)
    {
      hit->trim_t_right = 0;
    }

  hit->internal_alignmentlength = hit->nwalignmentlength
    - hit->trim_q_left - hit->trim_t_left
    - hit->trim_q_right - hit->trim_t_right;

  hit->internal_indels = hit->nwindels
    - hit->trim_q_left - hit->trim_t_left
    - hit->trim_q_right - hit->trim_t_right;

  hit->internal_gaps = hit->nwgaps
    - ((hit->trim_q_left  + hit->trim_t_left)  > 0 ? 1 : 0)
    - ((hit->trim_q_right + hit->trim_t_right) > 0 ? 1 : 0);

  /* CD-HIT */
  hit->id0 = hit->shortest > 0 ? 100.0 * hit->matches / hit->shortest : 0.0;
  /* all diffs */
  hit->id1 = hit->nwalignmentlength > 0 ?
    100.0 * hit->matches / hit->nwalignmentlength : 0.0;
  /* internal diffs */
  hit->id2 = hit->internal_alignmentlength > 0 ?
    100.0 * hit->matches / hit->internal_alignmentlength : 0.0;
  /* Marine Biology Lab */
  hit->id3 = std::max(0.0, 100.0 * (1.0 - (1.0 * (hit->mismatches + hit->nwgaps) /
                                      hit->longest)));
  /* BLAST */
  hit->id4 = hit->nwalignmentlength > 0 ?
    100.0 * hit->matches / hit->nwalignmentlength : 0.0;

  switch (opts->opt_iddef)
    {
    case 0:
      hit->id = hit->id0;
      break;
    case 1:
      hit->id = hit->id1;
      break;
    case 2:
      hit->id = hit->id2;
      break;
    case 3:
      hit->id = hit->id3;
      break;
    case 4:
      hit->id = hit->id4;
      break;
    }
}

/******************************************************************************/
auto Vsearch_Searchcore::acceptable_unaligned(struct searchinfo_s const & searchinfo,
                                 int const target,
                                 Vsearch_Database* db) -> bool
{
  /* consider whether a hit satisfies accepted criteria before alignment */

  // true: needs further consideration
  // false: reject

  auto const * qseq = searchinfo.qsequence;
  auto const * dlabel = db->getheader(target);
  auto const * dseq = db->getsequence(target);
  int64_t const dseqlen = db->getsequencelen(target);
  int64_t const tsize = db->getabundance(target);

  return (
          /* maxqsize */
          (searchinfo.qsize <= opts->opt_maxqsize)
          and
          /* mintsize */
          (tsize >= opts->opt_mintsize)
          and
          /* minsizeratio */
          (searchinfo.qsize >= opts->opt_minsizeratio * tsize)
          and
          /* maxsizeratio */
          (searchinfo.qsize <= opts->opt_maxsizeratio * tsize)
          and
          /* minqt */
          (searchinfo.qseqlen >= opts->opt_minqt * dseqlen)
          and
          /* maxqt */
          (searchinfo.qseqlen <= opts->opt_maxqt * dseqlen)
          and
          /* minsl */
          (searchinfo.qseqlen < dseqlen ?
           searchinfo.qseqlen >= opts->opt_minsl * dseqlen :
           dseqlen >= opts->opt_minsl * searchinfo.qseqlen)
          and
          /* maxsl */
          (searchinfo.qseqlen < dseqlen ?
           searchinfo.qseqlen <= opts->opt_maxsl * dseqlen :
           dseqlen <= opts->opt_maxsl * searchinfo.qseqlen)
          and
          /* idprefix */
          ((searchinfo.qseqlen >= opts->opt_idprefix) and
           (dseqlen >= opts->opt_idprefix) and
           (seqcmp.seqcmp(qseq, dseq, opts->opt_idprefix) == 0))
          and
          /* idsuffix */
          ((searchinfo.qseqlen >= opts->opt_idsuffix) and
           (dseqlen >= opts->opt_idsuffix) and
           (seqcmp.seqcmp(qseq + searchinfo.qseqlen - opts->opt_idsuffix,
                    dseq + dseqlen - opts->opt_idsuffix,
                    opts->opt_idsuffix) == 0))
          and
          /* self */
          ((opts->opt_self == 0) or (std::strcmp(searchinfo.query_head, dlabel) != 0))
          and
          /* selfid */
          ((opts->opt_selfid == 0) or
           (searchinfo.qseqlen != dseqlen) or
           (seqcmp.seqcmp(qseq, dseq, searchinfo.qseqlen) != 0))
          );
}

/******************************************************************************/
auto Vsearch_Searchcore::acceptable_aligned(struct searchinfo_s const & searchinfo,
                               struct hit * hit,
                               Vsearch_Database* db) -> bool
{
  if (/* weak_id */
      (hit->id >= 100.0 * opts->opt_weak_id) and
      /* maxsubs */
      (hit->mismatches <= opts->opt_maxsubs) and
      /* maxgaps */
      (hit->internal_gaps <= opts->opt_maxgaps) and
      /* mincols */
      (hit->internal_alignmentlength >= opts->opt_mincols) and
      /* leftjust */
      ((opts->opt_leftjust == 0) or (hit->trim_q_left +
                           hit->trim_t_left == 0)) and
      /* rightjust */
      ((opts->opt_rightjust == 0) or (hit->trim_q_right +
                            hit->trim_t_right == 0)) and
      /* query_cov */
      (hit->matches + hit->mismatches >= opts->opt_query_cov * searchinfo.qseqlen) and
      /* target_cov */
      (hit->matches + hit->mismatches >=
      opts->opt_target_cov * db->getsequencelen(hit->target)) and
      /* maxid */
      (hit->id <= 100.0 * opts->opt_maxid) and
      /* mid */
      (100.0 * hit->matches / (hit->matches + hit->mismatches) >= opts->opt_mid) and
      /* maxdiffs */
      (hit->mismatches + hit->internal_indels <= opts->opt_maxdiffs))
    {

      if (hit->id >= 100.0 * opts->opt_id)
        {
          /* accepted */
          hit->accepted = true;
          hit->weak = false;
          return true;
        }
      /* rejected, but weak hit */
      hit->rejected = true;
      hit->weak = true;
      return false;
    }

  /* rejected */
  hit->rejected = true;
  hit->weak = false;
  return false;
}

/******************************************************************************/
auto Vsearch_Searchcore::align_delayed(struct searchinfo_s * searchinfo,
                                       Vsearch_Database* db,
                                       Vsearch_AlignSimd* align) -> void
{
  /* compute global alignment */

  std::array<unsigned int, MAXDELAYED> target_list {{}};
  std::array<CELL, MAXDELAYED> nwscore_list {{}};
  std::array<unsigned short, MAXDELAYED> nwalignmentlength_list {{}};
  std::array<unsigned short, MAXDELAYED> nwmatches_list {{}};
  std::array<unsigned short, MAXDELAYED> nwmismatches_list {{}};
  std::array<unsigned short, MAXDELAYED> nwgaps_list {{}};
  std::array<char *, MAXDELAYED> nwcigar_list {{}};

  int target_count = 0;

  for (int x = searchinfo->finalized; x < searchinfo->hit_count; x++)
    {
      struct hit const * hit = searchinfo->hits + x;
      if (not hit->rejected)
        {
          target_list[target_count++] = hit->target;
        }
    }
  if (target_count != 0)
    {
      align->search16(searchinfo->s,
               target_count,
               target_list.data(),
               nwscore_list.data(),
               nwalignmentlength_list.data(),
               nwmatches_list.data(),
               nwmismatches_list.data(),
               nwgaps_list.data(),
               nwcigar_list.data(), db);
    }

  int i = 0;

  for (int x = searchinfo->finalized; x < searchinfo->hit_count; x++)
    {
      /* maxrejects or maxaccepts reached - ignore remaining hits */
      if ((searchinfo->rejects < opts->opt_maxrejects) and (searchinfo->accepts < opts->opt_maxaccepts))
        {
          struct hit * hit = searchinfo->hits + x;

          if (hit->rejected)
            {
              searchinfo->rejects++;
            }
          else
            {
              int64_t const target = hit->target;
              int64_t nwscore = nwscore_list[i];

              char * nwcigar = nullptr;
              int64_t nwalignmentlength = 0;
              int64_t nwmatches = 0;
              int64_t nwmismatches = 0;
              int64_t nwgaps = 0;

              int64_t const dseqlen = db->getsequencelen(target);

              if (nwscore == std::numeric_limits<short>::max())
                {
                  /* In case the SIMD aligner cannot align,
                     perform a new alignment with the
                     linear memory aligner */

                  char * dseq = db->getsequence(target);

                  if (nwcigar_list[i] != nullptr)
                    {
                      util.xfree(nwcigar_list[i]);
                    }

                  nwcigar = util.xstrdup(searchinfo->lma->align(searchinfo->qsequence,
                                                   dseq,
                                                   searchinfo->qseqlen,
                                                   dseqlen));

                  searchinfo->lma->alignstats(nwcigar,
                                      searchinfo->qsequence,
                                      dseq,
                                      & nwscore,
                                      & nwalignmentlength,
                                      & nwmatches,
                                      & nwmismatches,
                                      & nwgaps);
                }
              else
                {
                  nwalignmentlength = nwalignmentlength_list[i];
                  nwmatches = nwmatches_list[i];
                  nwmismatches = nwmismatches_list[i];
                  nwgaps = nwgaps_list[i];
                  nwcigar = nwcigar_list[i];
                }

              hit->aligned = true;
              hit->shortest = std::min(searchinfo->qseqlen, static_cast<int>(dseqlen));
              hit->longest = std::max(searchinfo->qseqlen, static_cast<int>(dseqlen));
              hit->nwalignment = nwcigar;
              hit->nwscore = nwscore;
              hit->nwdiff = nwalignmentlength - nwmatches;
              hit->nwgaps = nwgaps;
              hit->nwindels = nwalignmentlength - nwmatches - nwmismatches;
              hit->nwalignmentlength = nwalignmentlength;
              hit->nwid = 100.0 * (nwalignmentlength - hit->nwdiff) /
                nwalignmentlength;
              hit->matches = nwalignmentlength - hit->nwdiff;
              hit->mismatches = hit->nwdiff - hit->nwindels;

              /* trim alignment and compute numbers excluding terminal gaps */
              align_trim(hit);

              /* test accept/reject criteria after alignment */
              if (acceptable_aligned(*searchinfo, hit, db))
                {

                  searchinfo->accepts++;
                }
              else
                {
                  searchinfo->rejects++;
                }

              ++i;
            }
        }
    }

  /* free ignored alignments */
  while (i < target_count)
    {
      util.xfree(nwcigar_list[i++]);
    }

  searchinfo->finalized = searchinfo->hit_count;
}

/******************************************************************************/
auto Vsearch_Searchcore::onequery(struct searchinfo_s * searchinfo,
                                  int seqmask,
                                  Vsearch_Database* db,
                                  Vsearch_DBIndex* dbindex,
                                  Vsearch_AlignSimd* align,
                                  Vsearch_Minheap* minheap,
                                  Vsearch_Unique* unique) -> void
{

  searchinfo->hit_count = 0;

  align->search16_qprep(searchinfo->s, searchinfo->qsequence, searchinfo->qseqlen);

  struct Scoring scoring;
  scoring.match = opts->opt_match;
  scoring.mismatch = opts->opt_mismatch;
  scoring.gap_open_query_interior = opts->opt_gap_open_query_interior;
  scoring.gap_extension_query_interior = opts->opt_gap_extension_query_interior;
  scoring.gap_open_query_left = opts->opt_gap_open_query_left;
  scoring.gap_open_target_left = opts->opt_gap_open_target_left;
  scoring.gap_open_query_interior = opts->opt_gap_open_query_interior;
  scoring.gap_open_target_interior = opts->opt_gap_open_target_interior;
  scoring.gap_open_query_right = opts->opt_gap_open_query_right;
  scoring.gap_open_target_right = opts->opt_gap_open_target_right;
  scoring.gap_extension_query_left = opts->opt_gap_extension_query_left;
  scoring.gap_extension_target_left = opts->opt_gap_extension_target_left;
  scoring.gap_extension_query_interior = opts->opt_gap_extension_query_interior;
  scoring.gap_extension_target_interior = opts->opt_gap_extension_target_interior;
  scoring.gap_extension_query_right = opts->opt_gap_extension_query_right;
  scoring.gap_extension_target_right = opts->opt_gap_extension_target_right;

  searchinfo->lma = new LinearMemoryAligner(scoring);

  /* extract unique kmer samples from query*/
  unique->count(searchinfo->uh, 8,
               searchinfo->qseqlen, searchinfo->qsequence,
               &searchinfo->kmersamplecount, &searchinfo->kmersample, seqmask);

  /* find database sequences with the most kmer hits */
  topscores(searchinfo, db, dbindex, minheap);

  /* analyse targets with the highest number of kmer hits */
  searchinfo->accepts = 0;
  searchinfo->rejects = 0;
  searchinfo->finalized = 0;

  int delayed = 0;

  while ((searchinfo->finalized + delayed < opts->opt_maxaccepts + opts->opt_maxrejects - 1) and
         (searchinfo->rejects < opts->opt_maxrejects) and
         (searchinfo->accepts < opts->opt_maxaccepts) and
         (not minheap->isempty(searchinfo->m)))
    {
      elem_t const e = minheap->poplast(searchinfo->m);

      struct hit * hit = searchinfo->hits + searchinfo->hit_count;

      hit->target = e.seqno;
      hit->count = e.count;
      hit->strand = searchinfo->strand;
      hit->rejected = false;
      hit->accepted = false;
      hit->aligned = false;
      hit->weak = false;
      hit->nwalignment = nullptr;

      /* Test some accept/reject criteria before alignment */
      if (acceptable_unaligned(*searchinfo, e.seqno, db))
        {
          ++delayed;
        }
      else
        {
          hit->rejected = true;
        }

      searchinfo->hit_count++;

      if (delayed == MAXDELAYED)
        {
          align_delayed(searchinfo, db, align);
          delayed = 0;
        }
     }

  if (delayed > 0)
    {
      align_delayed(searchinfo, db, align);
    }

  delete searchinfo->lma;
}

/******************************************************************************/
auto Vsearch_Searchcore::joinhits(struct searchinfo_s * si_plus,
                     struct searchinfo_s * si_minus,
                     std::vector<struct hit> & hits) -> void
{
  /* join and sort accepted and weak hits from both strands */
  /* free the remaining alignments */

  auto const counter = count_number_of_hits_to_keep(si_plus) + count_number_of_hits_to_keep(si_minus);

  /* allocate new array of hits */
  hits.reserve(counter);

  copy_over_hits_to_be_kept(hits, si_plus);
  copy_over_hits_to_be_kept(hits, si_minus);

  free_rejected_alignments(si_plus);
  free_rejected_alignments(si_minus);

  /* last, sort the hits */
  std::qsort(hits.data(), counter, sizeof(struct hit), hit_compare_byid);
}
/******************************************************************************/
