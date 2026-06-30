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

#include "vsearch.h"
#include "bitmap.h"
#include "dbindex.h"

#include <array>
#include <cstdint>  // uint64_t
#include <cstdio>  // std::FILE, std::fprintf
#include <cstring>  // std::memset
#include <iterator>  // std::next

/******************************************************************************/
auto Vsearch_DBIndex::getbitmap(unsigned int const kmer) -> unsigned char *
{
  auto * a_bitmap_s = *std::next(kmerbitmap, kmer);
  if (a_bitmap_s != nullptr)
    {
      return a_bitmap_s->bitmap;
    }
  return nullptr;
}
/******************************************************************************/

auto Vsearch_DBIndex::getmatchcount(unsigned int const kmer) -> unsigned int
{
  return *std::next(kmercount, kmer);
}

/******************************************************************************/
auto Vsearch_DBIndex::getmatchlist(unsigned int const kmer) -> unsigned int *
{
  return std::next(kmerindex, *std::next(kmerhash, kmer));
}
/******************************************************************************/

auto Vsearch_DBIndex::getmapping(unsigned int const index) -> unsigned int
{
  return *std::next(dbindex_map, index);
}
/******************************************************************************/

auto Vsearch_DBIndex::getcount() -> unsigned int
{
  return dbindex_count;
}
/******************************************************************************/
auto Vsearch_DBIndex::addsequence(unsigned int seqno) -> void
{

  unsigned int uniquecount = 0;
  unsigned int const * uniquelist = nullptr;
  unique->count(dbindex_uh, 8,
               db->getsequencelen(seqno), db->getsequence(seqno),
               &uniquecount, &uniquelist, 1);
  dbindex_map[dbindex_count] = seqno;
  for (auto i = 0U; i < uniquecount; i++)
    {
      auto const kmer = uniquelist[i];
      if (kmerbitmap[kmer] != nullptr)
        {
          ++kmercount[kmer];
          bitmap->set(kmerbitmap[kmer], dbindex_count);
        }
      else
        {
          kmerindex[kmerhash[kmer] + kmercount[kmer]] = dbindex_count;
          ++kmercount[kmer];
        }
    }
  ++dbindex_count;
}

/******************************************************************************/
auto Vsearch_DBIndex::addallsequences() -> void
{
  unsigned int const seqcount = db->getsequencecount();
  for (auto seqno = 0U; seqno < seqcount ; seqno++)
    {
      addsequence(seqno);
    }
}

/******************************************************************************/
auto Vsearch_DBIndex::prepare() -> void
{
  dbindex_uh = unique->init();

  unsigned int const seqcount = db->getsequencecount();
  kmerhashsize = 1U << (2 * 8);

  /* allocate memory for kmer count array */
  kmercount = (unsigned int *) util.xmalloc(kmerhashsize * sizeof(unsigned int));
  std::memset(kmercount, 0, kmerhashsize * sizeof(unsigned int));

  /* first scan, just count occurences */
  //progress_init("Counting k-mers", seqcount);
  for (auto seqno = 0U; seqno < seqcount ; seqno++)
    {
      unsigned int uniquecount = 0;
      unsigned int const * uniquelist = nullptr;
      unique->count(dbindex_uh, 8,
                   db->getsequencelen(seqno), db->getsequence(seqno),
                   &uniquecount, &uniquelist, 1);
      for (auto i = 0U; i < uniquecount; i++)
        {
          ++kmercount[uniquelist[i]];
        }
    }

  // use_bitmap == 1
  bitmap_mincount = seqcount / bitmap_threshold;

  /* allocate and zero bitmap pointers */
  kmerbitmap = (struct bitmap_s **) util.xmalloc(kmerhashsize * sizeof(struct bitmap_s *));
  std::memset(kmerbitmap, 0, kmerhashsize * sizeof(struct bitmap_s *));

  /* hash / bitmap setup */
  /* convert hash counts to position in index */
  kmerhash = (uint64_t *) util.xmalloc((kmerhashsize + 1) * sizeof(uint64_t));
  uint64_t sum = 0;
  for (auto i = 0U; i < kmerhashsize; i++)
    {
      kmerhash[i] = sum;
      if (kmercount[i] >= bitmap_mincount)
        {
          // Calculate a size that is cleanly aligned to 16-byte blocks
          // to completely satisfy ASan container-overflow checks.
          // pad for xmm
          unsigned int aligned_seqcount = ((seqcount + 127 + 15) / 16) * 16;
          kmerbitmap[i] = bitmap->init(aligned_seqcount);
          bitmap->reset_all(kmerbitmap[i]);
        }
      else
        {
          sum += kmercount[i];
        }
    }
  kmerindexsize = sum;
  kmerhash[kmerhashsize] = sum;

  /* reset counts */
  std::memset(kmercount, 0, kmerhashsize * sizeof(unsigned int));

  /* allocate space for actual data */
  kmerindex = (unsigned int *) util.xmalloc(kmerindexsize * sizeof(unsigned int));

  /* allocate space for mapping from indexno to seqno */
  dbindex_map = (unsigned int *) util.xmalloc(seqcount * sizeof(unsigned int));

  dbindex_count = 0;
}
/******************************************************************************/
