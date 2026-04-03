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

#include "align_simd.h"
#include "vsearch_main.h"
#include "dbindex.h"
#include "linmemalign.h"
#include "minheap.h"
#include "unique.h"
#include "span.h"
#include <algorithm>  // std::copy, std::fill, std::fill_n, std::max, std::max_element, std::min, std::transform
#include <array>
#include <cassert>
#include <cctype>  // std::tolower
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint> // int64_t, uint64_t
#include <cstdlib>  // std::qsort
#include <cstdio>  // std::FILE, std::fprintf, std::sscanf
#include <cstring>  // std::strlen, std::strcpy
#include <iterator>  // std::next
#include <limits>
#include <numeric>  // std::accumulate
#include <pthread.h>
#include <vector>


/*
  This code implements the method described in this paper:

  Robert C. Edgar, Brian J. Haas, Jose C. Clemente, Christopher Quince
  and Rob Knight (2011)
  UCHIME improves sensitivity and speed of chimera detection
  Bioinformatics, 27, 16, 2194-2200
  https://doi.org/10.1093/bioinformatics/btr381
*/

/******************************************************************************/
Vsearch_Main::Vsearch_Main() {
    opts = Vsearch_Options::getInstance();
    seqno = 0;

    align = new Vsearch_AlignSimd();
    search = new Vsearch_Searchcore();
    minheap = new Vsearch_Minheap();
    unique = new Vsearch_Unique();
    bitmap = new Vsearch_Bitmap();
    cigar = new Vsearch_Cigar();
    maps = new Vsearch_Maps();

    chimera_count = 0;
    nonchimera_count = 0;
    borderline_count = 0;
    total_count = 0;
    chimera_abundance = 0;
    nonchimera_abundance = 0;
    borderline_abundance = 0;
    total_abundance = 0;
    parts = 0;

    cpu_features_detect();

    /* defaults */
    opts->opt_maxaccepts = few;
    opts->opt_maxrejects = rejects;
    opts->opt_id = chimera_id;
}
/******************************************************************************/
Vsearch_Main::~Vsearch_Main() {
    delete search;
    delete align;
    delete minheap;
    delete unique;
    delete bitmap;
    delete cigar;
    delete maps;
}
/******************************************************************************/
auto Vsearch_Main::cpu_features_detect() -> void {
#ifdef __aarch64__
#ifdef __ARM_NEON
    /* may check /proc/cpuinfo for asimd or neon */
    opts->neon_present = 1;
#else
#error ARM Neon not present
#endif
#elif __PPC__
    altivec_present = 1;
#elif __x86_64__
    unsigned int a = 0;
    unsigned int b = 0;
    unsigned int c = 0;
    unsigned int d = 0;

    cpuid(0, 0, a, b, c, d);
    unsigned int const maxlevel = a & 0xff;

    if (maxlevel >= 1)
    {
        cpuid(1, 0, a, b, c, d);
        //mmx_present    = (d >> 23U) & 1U;
        //sse_present    = (d >> 25U) & 1U;
        opts->sse2_present   = (d >> 26U) & 1U;
        //sse3_present   = (c >>  0U) & 1U;
        opts->ssse3_present  = (c >>  9U) & 1U;
        //sse41_present  = (c >> 19U) & 1U;
        //sse42_present  = (c >> 20U) & 1U;
        //popcnt_present = (c >> 23U) & 1U;
       // avx_present    = (c >> 28U) & 1U;

        if (maxlevel >= 7)
        {
            cpuid(7, 0, a, b, c, d);
            //avx2_present = (b >>  5U) & 1U;
        }
    }
#else
    // simde
#endif
}
/******************************************************************************/
auto Vsearch_Main::realloc_arrays(struct chimera_info_s * chimera_info) -> void
{

  /* default for uchime, uchime2, and uchime3 */
  parts = 4;

  const int maxhlen = std::max(chimera_info->query_head_len, 1);
  if (maxhlen > chimera_info->head_alloc)
    {
      chimera_info->query_head.resize(maxhlen + 1);
    }
  chimera_info->head_alloc = std::max(chimera_info->head_alloc, maxhlen);

  /* realloc arrays based on query length */

  const int maxqlen = std::max(chimera_info->query_len, 1);
  auto const max_2x2_size = static_cast<size_t>(maxcandidates * maxqlen);

  if (maxqlen > chimera_info->query_alloc)
    {
      chimera_info->query_alloc = maxqlen;

      chimera_info->query_seq.resize(maxqlen + 1);

      chimera_info->maxi.resize(maxqlen + 1);
      chimera_info->maxsmooth.resize(maxqlen);
      chimera_info->match.resize(max_2x2_size);
      chimera_info->insert.resize(max_2x2_size);
      chimera_info->smooth.resize(max_2x2_size);

      chimera_info->scan_p.resize(maxqlen + 1);
      chimera_info->scan_q.resize(maxqlen + 1);

      const int maxalnlen = maxqlen + (2 * db->getlongestsequence());
      chimera_info->paln.resize(maxparents);
      for (auto & a_parent_alignment : chimera_info->paln) {
        a_parent_alignment.resize(maxalnlen + 1);
      }
      chimera_info->qaln.resize(maxalnlen + 1);
      chimera_info->diffs.resize(maxalnlen + 1);
      chimera_info->votes.resize(maxalnlen + 1);
      chimera_info->model.resize(maxalnlen + 1);
      chimera_info->ignore.resize(maxalnlen + 1);
    }

  // resize query parts if longer than earlier, minimum 100
  const int maxpartlen = std::max((maxqlen + parts - 1) / parts, 100);
  if (maxpartlen > chimera_info->part_alloc)
    {
      for (auto & query_info: chimera_info->si)
        {
          query_info.qsequence_v.resize(maxpartlen + 1);
          query_info.qsequence = query_info.qsequence_v.data();
        }
      chimera_info->part_alloc = maxpartlen;
    }
}
/******************************************************************************/
auto Vsearch_Main::reset_matches(struct chimera_info_s * a_chimera_info) -> void {
  // refactoring: initialization to zero? (useless), or reset to zero??
  std::fill(a_chimera_info->match.begin(), a_chimera_info->match.end(), 0);
  std::fill(a_chimera_info->insert.begin(), a_chimera_info->insert.end(), 0);
}
/******************************************************************************/
auto Vsearch_Main::find_matches(struct chimera_info_s * chimera_info) -> void
{
  /* find the positions with matches for each potential parent */
  /* also note the positions with inserts in front */

  auto & qseq = chimera_info->query_seq;

  for (auto i = 0; i < chimera_info->cand_count; ++i)
    {
      auto const * tseq = db->getsequence(chimera_info->cand_list[i]);

      auto qpos = 0;
      auto tpos = 0;

      auto * cigar_start = chimera_info->nwcigar[i];
      auto const cigar_length = std::strlen(cigar_start);
      auto const cigar_pairs = cigar->parse_cigar_string(Span<char>{cigar_start, cigar_length});

      for (auto const & a_pair: cigar_pairs) {
        auto const operation = a_pair.first;
        auto const runlength = a_pair.second;
        switch (operation) {
        case Operation::match:
          for (auto j = 0; j < runlength; ++j)
            {
              if ((maps->map_4bit(qseq[qpos]) &
                  maps->map_4bit(tseq[tpos])) != 0U)
                {
                  chimera_info->match[(i * chimera_info->query_len) + qpos] = 1;
                }
              ++qpos;
              ++tpos;
            }
          break;

        case Operation::insertion:
          chimera_info->insert[(i * chimera_info->query_len) + qpos] = runlength;
          tpos += runlength;
          break;

        case Operation::deletion:
          qpos += runlength;
          break;
        }
      }
    }
}
/******************************************************************************/
auto Vsearch_Main::find_best_parents(struct chimera_info_s * ci) -> int
{
  reset_matches(ci);
  find_matches(ci);

  std::array<int, maxparents> best_parent_cand {{}};

  for (int f = 0; f < 2; ++f)
    {
      best_parent_cand[f] = -1;
      ci->best_parents[f] = -1;
    }

  std::vector<bool> cand_selected(ci->cand_count, false);

  for (int f = 0; f < 2; ++f)
    {
      if (f > 0)
        {
          /* for all parents except the first */

          /* wipe out matches for all candidates in positions
             covered by the previous parent */

          for (int qpos = window - 1; qpos < ci->query_len; ++qpos)
            {
              int const z = (best_parent_cand[f - 1] * ci->query_len) + qpos;
              if (ci->smooth[z] == ci->maxsmooth[qpos])
                {
                  for (int i = qpos + 1 - window; i <= qpos; ++i)
                    {
                      for (int j = 0; j < ci->cand_count; ++j)
                        {
                          ci->match[(j * ci->query_len) + i] = 0;
                        }
                    }
                }
            }
        }


      /* Compute smoothed score in a 32bp window for each candidate. */
      /* Record max smoothed score for each position among candidates left. */

      // refactoring: reset or initialize?
      std::fill(ci->maxsmooth.begin(), ci->maxsmooth.end(), 0);

      for (auto i = 0; i < ci->cand_count; ++i)
        {
          if (not cand_selected[i])
            {
              int sum = 0;
              for (int qpos = 0; qpos < ci->query_len; ++qpos)
                {
                  int const z = (i * ci->query_len) + qpos;
                  sum += ci->match[z];

                  if (qpos >= window)
                    {
                      sum -= ci->match[z - window];
                    }
                  if (qpos >= window - 1)
                    {
                      ci->smooth[z] = sum;
                      ci->maxsmooth[qpos] = std::max(ci->smooth[z], ci->maxsmooth[qpos]);
                    }
                }
            }
        }

      /* find parent with the most wins */

      std::vector<int> wins(ci->cand_count, 0);

      for (int qpos = window - 1; qpos < ci->query_len; ++qpos)
        {
          if (ci->maxsmooth[qpos] != 0)
            {
              for (auto i = 0; i < ci->cand_count; ++i)
                {
                  if (not cand_selected[i])
                    {
                      int const z = (i * ci->query_len) + qpos;
                      if (ci->smooth[z] == ci->maxsmooth[qpos])
                        {
                          ++wins[i];
                        }
                    }
                }
            }
        }

      /* select best parent based on most wins */

      int maxwins = 0;
      for (auto i = 0; i < ci->cand_count; ++i)
        {
          int const w = wins[i];
          if (w > maxwins)
            {
              maxwins = w;
              best_parent_cand[f] = i;
            }
        }

      /* terminate loop if no parent found */

      if (best_parent_cand[f] < 0) {
        break;
      }

#if 0
      printf("Query %d: Best parent (%d) candidate: %d. Wins: %d\n",
             ci->query_no, f, best_parent_cand[f], maxwins);
#endif

      ci->best_parents[f] = best_parent_cand[f];
      cand_selected[best_parent_cand[f]] = true;
    }

  /* Check if at least 2 candidates selected */

  return static_cast<int>((best_parent_cand[0] >= 0) and (best_parent_cand[1] >= 0));
}
/******************************************************************************/
auto Vsearch_Main::find_total_alignment_length(struct chimera_info_s const * chimera_info) -> int
{
  // query_len, plus the sum of the longest insertion runs (I) for each position
  return std::accumulate(chimera_info->maxi.begin(),
                         chimera_info->maxi.end(),
                         chimera_info->query_len);
}
/******************************************************************************/
auto Vsearch_Main::fill_max_alignment_length(struct chimera_info_s * chimera_info) -> void
{
  /* find max insertions in front of each position in the query sequence */

  std::fill(chimera_info->maxi.begin(), chimera_info->maxi.end(), 0);

  auto const count = static_cast<size_t>(chimera_info->parents_found);
  assert(count <= chimera_info->best_parents.size());
  auto const best_parents_view = Span<int>{chimera_info->best_parents.data(), count};
  for (auto const best_parent : best_parents_view) {
    auto pos = 0LL;
    auto * cigar_start = chimera_info->nwcigar[best_parent];
    auto const cigar_length = std::strlen(cigar_start);
    auto const cigar_pairs = cigar->parse_cigar_string(Span<char>{cigar_start, cigar_length});

    for (auto const & a_pair: cigar_pairs) {
      auto const operation = a_pair.first;
      auto const runlength = a_pair.second;
      switch (operation) {
      case Operation::match:
      case Operation::deletion:
        pos += runlength;
        break;

      case Operation::insertion:
        assert(runlength <= std::numeric_limits<int>::max());
        chimera_info->maxi[pos] = std::max(static_cast<int>(runlength), chimera_info->maxi[pos]);
        break;
      }
    }
  }
}
/******************************************************************************/
auto Vsearch_Main::fill_alignment_parents(struct chimera_info_s * ci) -> void
{
  /* fill in alignment strings for the parents */

  for (auto i = 0; i < ci->parents_found; ++i)
    {
      auto & alignment = ci->paln[i];
      int const cand = ci->best_parents[i];
      int const target_seqno = ci->cand_list[cand];
      char const * target_seq = db->getsequence(target_seqno);

      auto is_inserted = false;
      int qpos = 0;
      int tpos = 0;
      int alnpos = 0;

      auto * cigar_start = ci->nwcigar[cand];
      auto const cigar_length = std::strlen(cigar_start);
      auto const cigar_pairs = cigar->parse_cigar_string(Span<char>{cigar_start, cigar_length});
      for (auto const & a_pair: cigar_pairs) {
        auto const operation = a_pair.first;
        auto const runlength = a_pair.second;
        switch (operation) {
        case Operation::insertion:
          for (int j = 0; j < ci->maxi[qpos]; ++j)
            {
              if (j < runlength)
                {
                  alignment[alnpos] = maps->map_uppercase(target_seq[tpos]);
                  ++tpos;
                  ++alnpos;
                }
              else
                {
                  alignment[alnpos] = '-';
                  ++alnpos;
                }
            }
          is_inserted = true;
          break;

        case Operation::match:
        case Operation::deletion:
          for (int j = 0; j < runlength; ++j)
            {
              if (not is_inserted)
                {
                  std::fill_n(&alignment[alnpos], ci->maxi[qpos], '-');
                  alnpos += ci->maxi[qpos];
                }

              if (operation == Operation::match)
                {
                  alignment[alnpos] = maps->map_uppercase(target_seq[tpos]);
                  ++tpos;
                  ++alnpos;
                }
              else
                {
                  alignment[alnpos] = '-';
                  ++alnpos;
                }

              ++qpos;
              is_inserted = false;
            }
        }
      }

      /* add any gaps at the end */

      if (not is_inserted)
        {
          std::fill_n(&alignment[alnpos], ci->maxi[qpos], '-');
          alnpos += ci->maxi[qpos];
        }

      /* end of sequence string */
      alignment[alnpos] = '\0';
    }
}
/******************************************************************************/
auto Vsearch_Main::eval_parents(struct chimera_info_s * ci) -> Status
{
  auto status = Status::no_alignment;
  ci->parents_found = 2;

  fill_max_alignment_length(ci);
  auto const alnlen = find_total_alignment_length(ci);

  fill_alignment_parents(ci);

  /* fill in alignment string for query */

  char * q = ci->qaln.data();
  int qpos = 0;
  for (auto i = 0; i < ci->query_len; ++i)
    {
      for (int j = 0; j < ci->maxi[i]; ++j)
        {
          *q = '-';
          ++q;
        }
      *q = maps->map_uppercase(ci->query_seq[qpos]);
      ++qpos;
      ++q;
    }
  for (int j = 0; j < ci->maxi[ci->query_len]; ++j)
    {
      *q = '-';
      ++q;
    }
  *q = 0;

  /* mark positions to ignore in voting */
  std::fill(ci->ignore.begin(), ci->ignore.end(), false);

  for (auto i = 0; i < alnlen; ++i)
    {
      auto const qsym  = maps->map_4bit(ci->qaln[i]);
      auto const p1sym = maps->map_4bit(ci->paln[0][i]);
      auto const p2sym = maps->map_4bit(ci->paln[1][i]);

      /* ignore gap positions and those next to the gap */
      if ((qsym == 0U) or (p1sym == 0U) or (p2sym == 0U))
        {
          ci->ignore[i] = true;
          if (i > 0)
            {
              ci->ignore[i - 1] = true;
            }
          if (i < alnlen - 1)
            {
              ci->ignore[i + 1] = true;
            }
        }

      /* ignore ambiguous symbols */
      if (maps->is_ambiguous_4bit(qsym) or
              maps->is_ambiguous_4bit(p1sym) or
              maps->is_ambiguous_4bit(p2sym))
        {
          ci->ignore[i] = true;
        }

      /* lower case parent symbols that differ from query */

      if ((p1sym != 0U) and (p1sym != qsym))
        {
          ci->paln[0][i] = std::tolower(ci->paln[0][i]);
        }

      if ((p2sym != 0U) and (p2sym != qsym))
        {
          ci->paln[1][i] = std::tolower(ci->paln[1][i]);
        }

      /* compute diffs */

      char diff = '\0';

      if ((qsym != 0U) and (p1sym != 0U) and (p2sym != 0U))
        {
          if (p1sym == p2sym)
            {
              if (qsym == p1sym)
                {
                  diff = ' ';
                }
              else
                {
                  diff = 'N';
                }
            }
          else
            {
              if (qsym == p1sym)
                {
                  diff = 'A';
                }
              else if (qsym == p2sym)
                {
                  diff = 'B';
                }
              else
                {
                  diff = '?';
                }
            }
        }
      else
        {
          diff = ' ';
        }

      ci->diffs[i] = diff;
    }

  ci->diffs[alnlen] = '\0';

  /* compute score */

  int sumA = 0;
  int sumB = 0;
  int sumN = 0;

  for (auto i = 0; i < alnlen; ++i)
    {
      if (ci->ignore[i]) { continue; }
      auto const diff = ci->diffs[i];

      if (diff == 'A')
        {
          ++sumA;
        }
      else if (diff == 'B')
        {
          ++sumB;
        }
      else if (diff != ' ')
        {
          ++sumN;
        }
    }

  int left_n = 0;
  int left_a = 0;
  int left_y = 0;
  int right_n = sumA;
  int right_a = sumN;
  int right_y = sumB;

  double best_h = -1;
  int best_i = -1;
  auto best_is_reverse = false;

  int best_left_y = 0;
  int best_right_y = 0;
  int best_left_n = 0;
  int best_right_n = 0;
  int best_left_a = 0;
  int best_right_a = 0;

  for (auto i = 0; i < alnlen; ++i)
    {
      if (not ci->ignore[i])
        {
          char const diff = ci->diffs[i];
          if (diff != ' ')
            {
              if (diff == 'A')
                {
                  ++left_y;
                  --right_n;
                }
              else if (diff == 'B')
                {
                  ++left_n;
                  --right_y;
                }
              else
                {
                  ++left_a;
                  --right_a;
                }

              double left_h = 0;
              double right_h = 0;
              double h = 0;

              if ((left_y > left_n) and (right_y > right_n))
                {
                  left_h = left_y / (opts->opt_xn * (left_n + opts->opt_dn) + left_a);
                  right_h = right_y / (opts->opt_xn * (right_n + opts->opt_dn) + right_a);
                  h = left_h * right_h;

                  if (h > best_h)
                    {
                      best_is_reverse = false;
                      best_h = h;
                      best_i = i;
                      best_left_n = left_n;
                      best_left_y = left_y;
                      best_left_a = left_a;
                      best_right_n = right_n;
                      best_right_y = right_y;
                      best_right_a = right_a;
                    }
                }
              else if ((left_n > left_y) and (right_n > right_y))
                {
                  /* swap left/right and yes/no */

                  left_h = left_n / (opts->opt_xn * (left_y + opts->opt_dn) + left_a);
                  right_h = right_n / (opts->opt_xn * (right_y + opts->opt_dn) + right_a);
                  h = left_h * right_h;

                  if (h > best_h)
                    {
                      best_is_reverse = true;
                      best_h = h;
                      best_i = i;
                      best_left_n = left_y;
                      best_left_y = left_n;
                      best_left_a = left_a;
                      best_right_n = right_y;
                      best_right_y = right_n;
                      best_right_a = right_a;
                    }
                }
            }
        }
    }

  ci->best_h = best_h > 0 ? best_h : 0.0;

  if (best_h >= 0.0)
    {
      status = Status::low_score;

      /* flip A and B if necessary */

      if (best_is_reverse)
        {
          for (auto i = 0; i < alnlen; ++i)
            {
              char const diff = ci->diffs[i];
              if (diff == 'A')
                {
                  ci->diffs[i] = 'B';
                }
              else if (diff == 'B')
                {
                  ci->diffs[i] = 'A';
                }
            }
        }

      /* fill in votes and model */

      for (auto i = 0; i < alnlen; ++i)
        {
          char const m = i <= best_i ? 'A' : 'B';
          ci->model[i] = m;

          char v = ' ';
          if (not ci->ignore[i])
            {
              char const d = ci->diffs[i];

              if ((d == 'A') or (d == 'B'))
                {
                  if (d == m)
                    {
                      v = '+';
                    }
                  else
                    {
                      v = '!';
                    }
                }
              else if ((d == 'N') or (d == '?'))
                {
                  v = '0';
                }
            }
          ci->votes[i] = v;

          /* lower case diffs for no votes */
          if (v == '!')
            {
              ci->diffs[i] = std::tolower(ci->diffs[i]);
            }
        }

      /* fill in crossover region */

      for (int i = best_i + 1; i < alnlen; ++i)
        {
          if ((ci->diffs[i] == ' ') or (ci->diffs[i] == 'A'))
            {
              ci->model[i] = 'x';
            }
          else
            {
              break;
            }
        }

      ci->votes[alnlen] = 0;
      ci->model[alnlen] = 0;

      /* count matches */

      auto const index_a = best_is_reverse ? 1U : 0U;
      auto const index_b = best_is_reverse ? 0U : 1U;

      int match_QA = 0;
      int match_QB = 0;
      int match_AB = 0;
      int match_QM = 0;
      int cols = 0;

      for (auto i = 0; i < alnlen; i++)
        {
          if (not ci->ignore[i])
            {
              ++cols;

              auto const qsym = maps->map_4bit(ci->qaln[i]);
              auto const asym = maps->map_4bit(ci->paln[index_a][i]);
              auto const bsym = maps->map_4bit(ci->paln[index_b][i]);
              auto const msym = (i <= best_i) ? asym : bsym;

              if (qsym == asym)
                {
                  ++match_QA;
                }

              if (qsym == bsym)
                {
                  ++match_QB;
                }

              if (asym == bsym)
                {
                  ++match_AB;
                }

              if (qsym == msym)
                {
                  ++match_QM;
                }
            }
        }

      int const seqno_a = ci->cand_list[ci->best_parents[index_a]];
      int const seqno_b = ci->cand_list[ci->best_parents[index_b]];

      double const QA = 100.0 * match_QA / cols;
      double const QB = 100.0 * match_QB / cols;
      double const AB = 100.0 * match_AB / cols;
      double const QT = std::max(QA, QB);
      double const QM = 100.0 * match_QM / cols;
      double const divdiff = QM - QT;

      int const sumL = best_left_n + best_left_a + best_left_y;
      int const sumR = best_right_n + best_right_a + best_right_y;

        if (best_h >= opts->opt_minh)
          {
            status = Status::suspicious;
            if ((divdiff >= opts->opt_mindiv) and
                (sumL >= opts->opt_mindiffs) and
                (sumR >= opts->opt_mindiffs))
              {
                status = Status::chimeric;
              }
          }

        ChimeHit2 thisResult;
        thisResult.Score = best_h;
        thisResult.QLabel = ci->query_head.data();
        thisResult.ALabel = db->getheader(seqno_a);
        thisResult.BLabel = db->getheader(seqno_b);

        thisResult.QA = QA;
        thisResult.QB = QB;
        thisResult.QC = AB;
        thisResult.QM = QM;
        thisResult.QT = QT;

        if (QA >= QB) {
            thisResult.CLabel = db->getheader(seqno_a);
        }else {
            thisResult.CLabel = db->getheader(seqno_b);
        }
        thisResult.LA = best_left_a;
        thisResult.LY = best_left_y;
        thisResult.LN = best_left_n;
        thisResult.RA = best_right_a;
        thisResult.RY = best_right_y;
        thisResult.RN = best_right_n;
        thisResult.Div = divdiff;
        thisResult.status = status == Status::chimeric ? 'Y' : (status == Status::low_score ? 'N' : '?');
        results.push_back(thisResult);
    }
  return status;
}
/******************************************************************************/
auto Vsearch_Main::query_init(struct searchinfo_s * search_info) -> void
{
  static constexpr auto overflow_padding = 16U;  // 16 * sizeof(short) = 32 bytes
  search_info->hits_v.resize(tophits);
  search_info->hits = search_info->hits_v.data();
  search_info->kmers_v.reserve(db->getsequencecount() + overflow_padding);
  search_info->kmers_v.resize(db->getsequencecount());
  search_info->kmers = search_info->kmers_v.data();
  search_info->hit_count = 0;
  search_info->uh = unique->init();
  search_info->s = align->search16_init(opts->opt_match,
                                 opts->opt_mismatch,
                                 opts->opt_gap_open_query_left,
                                 opts->opt_gap_open_target_left,
                                 opts->opt_gap_open_query_interior,
                                 opts->opt_gap_open_target_interior,
                                 opts->opt_gap_open_query_right,
                                 opts->opt_gap_open_target_right,
                                 opts->opt_gap_extension_query_left,
                                 opts->opt_gap_extension_target_left,
                                 opts->opt_gap_extension_query_interior,
                                 opts->opt_gap_extension_target_interior,
                                 opts->opt_gap_extension_query_right,
                                 opts->opt_gap_extension_target_right,
                                 opts);
  search_info->m = minheap->init(tophits);
}

/******************************************************************************/
auto Vsearch_Main::query_exit(struct searchinfo_s * search_info) -> void
{
  align->search16_exit(search_info->s);
  unique->exit(search_info->uh);
  minheap->exit(search_info->m);

  search_info->qsequence = nullptr;
  search_info->hits = nullptr;
  search_info->kmers = nullptr;
}
/******************************************************************************/
auto Vsearch_Main::partition_query(struct chimera_info_s * chimera_info) -> void
{
  auto rest = chimera_info->query_len;
  auto * cursor = chimera_info->query_seq.data();
  for (auto i = 0; i < parts; ++i)
    {
      auto const length = (rest + (parts - i - 1)) / (parts - i);

      auto & search_info = chimera_info->si[i];

      search_info.query_no = chimera_info->query_no;
      search_info.strand = 0;
      search_info.qsize = chimera_info->query_size;
      search_info.query_head_len = chimera_info->query_head_len;
      search_info.query_head = chimera_info->query_head.data();
      search_info.qseqlen = length;
      assert(static_cast<std::size_t>(length) <= search_info.qsequence_v.size());
      std::copy(cursor, std::next(cursor, length), search_info.qsequence_v.begin());
      search_info.qsequence_v[length] = '\0';

      rest -= length;
      cursor = std::next(cursor, length);
    }
}
/******************************************************************************/
auto Vsearch_Main::chimera_init(struct chimera_info_s * ci) -> void
{

  for (auto i = 0; i < maxparts; ++i)
    {
      query_init(&ci->si[i]);
    }

  ci->s = align->search16_init(opts->opt_match,
                        opts->opt_mismatch,
                        opts->opt_gap_open_query_left,
                        opts->opt_gap_open_target_left,
                        opts->opt_gap_open_query_interior,
                        opts->opt_gap_open_target_interior,
                        opts->opt_gap_open_query_right,
                        opts->opt_gap_open_target_right,
                        opts->opt_gap_extension_query_left,
                        opts->opt_gap_extension_target_left,
                        opts->opt_gap_extension_query_interior,
                        opts->opt_gap_extension_target_interior,
                        opts->opt_gap_extension_query_right,
                        opts->opt_gap_extension_target_right, opts);
}
/******************************************************************************/
auto Vsearch_Main::chimera_exit(struct chimera_info_s * ci) -> void
{
    align->search16_exit(ci->s);

  for (auto & a_search_info : ci->si) {
    query_exit(&a_search_info);
  }
}
/******************************************************************************/
auto Vsearch_Main::chimera_core(struct chimera_info_s * ci) -> uint64_t
{
  chimera_init(ci);

  std::vector<struct hit> allhits_list(maxcandidates);

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

  LinearMemoryAligner lma(scoring);

  while (true)
    {
      /* get next sequence */
      if (!uchimeDeNovo)
        {
          if (seqno < referenceThreadStop)
            {
              ci->query_head_len = (*queryNames)[seqno].length();
              ci->query_len = (*querySeqs)[seqno].length();
              ci->query_no = seqno;
              ci->query_size = (*queryAbunds)[seqno];

              /* if necessary expand memory for arrays based on query length */
              realloc_arrays(ci);

              /* copy the data locally (query seq, head) */
              std::strcpy(ci->query_head.data(), (*queryNames)[seqno].c_str());
              std::strcpy(ci->query_seq.data(), (*querySeqs)[seqno].c_str());
            }
          else
            {
              break; /* end while loop */
            }
        }
      else
        {
          if (seqno < db->getsequencecount())
            {
              ci->query_no = seqno;
              ci->query_head_len = db->getheaderlen(seqno);
              ci->query_len = db->getsequencelen(seqno);
              ci->query_size = db->getabundance(seqno);

              /* if necessary expand memory for arrays based on query length */
              realloc_arrays(ci);

              std::strcpy(ci->query_head.data(), db->getheader(seqno));
              std::strcpy(ci->query_seq.data(), db->getsequence(seqno));
            }
          else
            {
              break; /* end while loop */
            }
        }

      auto status = Status::no_parents;

      /* partition query */
      partition_query(ci);

      /* perform searches and collect candidate parents */
      ci->cand_count = 0;
      ci->best_h = 0.0;
      auto allhits_count = 0;

      if (ci->query_len >= parts)
        {
          std::vector<struct hit> hits;
          for (auto i = 0; i < parts; ++i)
            {
              search->onequery(&ci->si[i], 1,
                               db, dbindex, align, minheap, unique);
              search->joinhits(&ci->si[i], nullptr, hits);

              for (auto & hit : hits) {
                if (hit.accepted)
                  {
                    allhits_list[allhits_count] = hit;
                    ++allhits_count;
                  }
                else
                  {
                    // Unallocate alignments for weak hits
                    if (hit.nwalignment)
                      {
                        util.xfree(hit.nwalignment);
                        hit.nwalignment = nullptr;
                      }
                  }
              }
              hits.clear();
            }
        }

      for (auto i = 0; i < allhits_count; ++i)
        {
          unsigned int const target = allhits_list[i].target;

          /* skip duplicates */
          auto k = 0;
          for (k = 0; k < ci->cand_count; ++k)
            {
              if (ci->cand_list[k] == target)
                {
                  break;
                }
            }

          if (k == ci->cand_count)
            {
              ci->cand_list[ci->cand_count] = target;
              ++ci->cand_count;
            }

          /* deallocate cigar */
          if (allhits_list[i].nwalignment != nullptr)
            {
              util.xfree(allhits_list[i].nwalignment);
              allhits_list[i].nwalignment = nullptr;
            }
        }


      /* align full query to each candidate */
      align->search16_qprep(ci->s, ci->query_seq.data(), ci->query_len);

      align->search16(ci->s,
               ci->cand_count,
               ci->cand_list.data(),
               ci->snwscore.data(),
               ci->snwalignmentlength.data(),
               ci->snwmatches.data(),
               ci->snwmismatches.data(),
               ci->snwgaps.data(),
               ci->nwcigar.data(),
               db);

      for (auto i = 0; i < ci->cand_count; ++i)
        {
          int64_t const target = ci->cand_list[i];
          int64_t nwscore = ci->snwscore[i];
          char * nwcigar = nullptr;
          int64_t nwalignmentlength = 0;
          int64_t nwmatches = 0;
          int64_t nwmismatches = 0;
          int64_t nwgaps = 0;

          if (nwscore == std::numeric_limits<short>::max())
            {
              /* In case the SIMD aligner cannot align,
                 perform a new alignment with the
                 linear memory aligner */

              auto * tseq = db->getsequence(target);
              int64_t const tseqlen = db->getsequencelen(target);

              if (ci->nwcigar[i] != nullptr)
                {
                  util.xfree(ci->nwcigar[i]);
                }

              nwcigar = util.xstrdup(lma.align(ci->query_seq.data(),
                                          tseq,
                                          ci->query_len,
                                          tseqlen));
              lma.alignstats(nwcigar,
                             ci->query_seq.data(),
                             tseq,
                             & nwscore,
                             & nwalignmentlength,
                             & nwmatches,
                             & nwmismatches,
                             & nwgaps);

              ci->nwcigar[i] = nwcigar;
              ci->nwscore[i] = nwscore;
              ci->nwalignmentlength[i] = nwalignmentlength;
              ci->nwmatches[i] = nwmatches;
              ci->nwmismatches[i] = nwmismatches;
              ci->nwgaps[i] = nwgaps;
            }
          else
            {
              ci->nwscore[i] = ci->snwscore[i];
              ci->nwalignmentlength[i] = ci->snwalignmentlength[i];
              ci->nwmatches[i] = ci->snwmatches[i];
              ci->nwmismatches[i] = ci->snwmismatches[i];
              ci->nwgaps[i] = ci->snwgaps[i];
            }
        }

      /* find the best pair of parents, then compute score for them */
          if (find_best_parents(ci) != 0)
            {
              status = eval_parents(ci);
            }
          else
            {
              status = Status::no_parents;
            }

      /* output results */

      ++total_count;
      total_abundance += ci->query_size;

      if (status == Status::chimeric)
        {
          ++chimera_count;
          chimera_abundance += ci->query_size;
          chimeraNames.insert(ci->query_head.data());
        }

      if (status == Status::suspicious)
        {
          ++borderline_count;
          borderline_abundance += ci->query_size;
        }

      if (status < Status::suspicious)
        {
          ++nonchimera_count;
          nonchimera_abundance += ci->query_size;

          /* output no parents, no chimeras */
          if (status < Status::low_score)
            {
              ChimeHit2 thisResult;
              thisResult.H = ci->best_h;
              thisResult.QLabel = ci->query_head.data();
              thisResult.status = "N";
              results.push_back(thisResult);
            }
        }

      if (status < Status::suspicious)
        {
          /* uchime_denovo: add non-chimeras to db */
          if (uchimeDeNovo) {
              dbindex->addsequence(seqno);
          }
        }

      for (auto i = 0; i < ci->cand_count; ++i)
        {
          if (ci->nwcigar[i] != nullptr)
            {
              util.xfree(ci->nwcigar[i]);
            }
        }

      ++seqno;
    }

  chimera_exit(ci);

  return 0;
}
/******************************************************************************/
auto Vsearch_Main::wo(int len, const char *s, int *beg, int *end) -> int
{
    static constexpr auto dust_word = 3;
    static constexpr auto word_count = 1U << (2U * dust_word);  // 64
    static constexpr auto bitmask = word_count - 1;
    const auto l1 = len - dust_word + 1 - 5; /* smallest possible region is 8 */
if (l1 < 0)
{
    return 0;
}

auto bestv = 0;
auto besti = 0;
auto bestj = 0;
std::array<int, word_count> counts {{}};
std::array<int, dust_window> words {{}};
auto word = 0;

for (auto j = 0; j < len; j++)
{
    word <<= 2U;
    word |= maps->map_2bit(s[j]);
    words[j] = word & bitmask;
}

for (auto i = 0; i < l1; i++)
{
    counts.fill(0);  // reset counts to zero

    auto sum = 0;

    for (auto j = dust_word - 1; j < len - i; j++)
    {
        word = words[i + j];
        const auto c = counts[word];
        if (c != 0)
        {
            sum += c;
            const auto v = 10 * sum / j;

            if (v > bestv)
            {
                bestv = v;
                besti = i;
                bestj = j;
            }
        }
        ++counts[word];
    }
}

*beg = besti;
*end = besti + bestj;

return bestv;
}
/******************************************************************************/
auto Vsearch_Main::dust(char * seq, int len) -> void
{
    static constexpr auto dust_level = 20;
    static constexpr auto half_dust_window = dust_window / 2;
    auto a = 0;
    auto b = 0;

    /* make a local copy of the original sequence */
    std::vector<char> local_seq(len + 1);
    strcpy(local_seq.data(), seq);

    for (auto i = 0; i < len; i += half_dust_window)
    {
        const auto l = (len > i + dust_window) ? dust_window : len - i;
        const auto v = wo(l, &local_seq[i], &a, &b);

        if (v > dust_level)
        {

            for (auto j = a + i; j <= b + i; j++)
            {
                seq[j] = local_seq[j] | 32U;  // check_5th_bit (0x20)
            }

            if (b < half_dust_window)
            {
                i += half_dust_window - b;
            }
        }
    }
}
/******************************************************************************/
auto Vsearch_Main::dust_all(uint64_t seqcount) -> void {
    for (auto i = 0; i < seqcount; i++) {
        dust(db->getsequence(i), db->getsequencelen(i));
    }
}
/******************************************************************************/
auto Vsearch_Main::vmain(std::vector<std::string>& sequenceNames,
                                    std::vector<std::string>& sequences,
                                    std::vector<std::string>& refNames,
                                    std::vector<std::string>& refSeqs,
                                    std::vector<float>& abunds,
                                    unsigned int start, unsigned int stop) -> std::vector<ChimeHit2> {

#ifdef __x86_64__
   if (opts->sse2_present == 0) {
      throw Rcpp::exception("Sorry, this program requires a cpu with SSE2.");
   }
#endif

  // are we running with a reference, or denovo
  uchimeDeNovo = (refNames.size() == 0);

  queryNames = &sequenceNames;
  querySeqs = &sequences;

  std::vector<int64_t> int64s;
  int64s.reserve(abunds.size());
  std::transform(abunds.begin(), abunds.end(),
                 std::back_inserter(int64s),
                 [](float f) { return static_cast<int64_t>(f); });
  queryAbunds = &int64s;

  if (uchimeDeNovo) {
      opts->opt_self = 1;
      opts->opt_selfid = 1;
      opts->opt_threads = 1;
      opts->opt_maxsizeratio = 1.0 / opts->opt_abskew;
  }

  tophits = opts->opt_maxaccepts + opts->opt_maxrejects;
  chimera_count = 0;
  nonchimera_count = 0;
  seqno = 0;
  referenceThreadStop = 0;

  struct chimera_info_s cia_single;
  cia = &cia_single;

  db = new Vsearch_Database();
  dbindex = new Vsearch_DBIndex(db, bitmap, unique);

  /* prepare queries / database */
  if (!uchimeDeNovo) {

      seqno = start;
      referenceThreadStop = stop;

      std::vector<std::int64_t> refAbunds(refSeqs.size(), 1);
      db->read(refNames, refSeqs, refAbunds);

      dust_all(db->getsequencecount());

      dbindex->prepare();
      dbindex->addallsequences();
  }
  else {
      db->read(sequenceNames, sequences, int64s);

      dust_all(db->getsequencecount());
      dbindex->prepare();
  }

  chimera_core(cia);

  delete db;
  delete dbindex;

  return results;
}
/******************************************************************************/

