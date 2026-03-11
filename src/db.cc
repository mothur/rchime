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
#include "db.h"
#include "maps.h"
#include <algorithm>  // std::min, std::max
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::fprintf, std::size_t
#include <cstdlib>  // std::qsort
#include <cstring>  // std::memcpy, std::strcmp
#include <limits>


constexpr uint64_t memchunk = 16777216;  // 2^24

/******************************************************************************/
auto Vsearch_Database::setinfo(bool new_is_fastq,
                uint64_t new_sequences,
                uint64_t new_nucleotides,
                uint64_t new_longest,
                uint64_t new_shortest,
                uint64_t new_longestheader) -> void
{

  sequences = new_sequences;
  nucleotides = new_nucleotides;
  longest = new_longest;
  shortest = new_shortest;
  longestheader = new_longestheader;
}
/******************************************************************************/
// private
auto Vsearch_Database::add(char const * header,
            char const * sequence,
            char const * quality,
            size_t const headerlength,
            size_t const sequencelength,
            int64_t const abundance) -> void
{
  /* Add a sequence to the database. Assumes that the database has been initialized. */

  /* grow space for data, if necessary */

  size_t const dataalloc_old = dataalloc;

  size_t needed = datalen + headerlength + 1 + sequencelength + 1;

  while (dataalloc < needed)
    {
      dataalloc += memchunk;
    }
  if (dataalloc > dataalloc_old)
    {
      datap = (char *) util.xrealloc(datap, dataalloc);
    }

  /* store the header */
  size_t const header_p = datalen;
  std::memcpy(datap + header_p,
              header,
              headerlength + 1);
  datalen += headerlength + 1;

  /* store sequence */
  size_t const sequence_p = datalen;
  std::memcpy(datap + sequence_p,
              sequence,
              sequencelength + 1);
  datalen += sequencelength + 1;

  size_t const quality_p = datalen;

  /* grow space for index, if necessary */
  size_t const seqindex_alloc_old = seqindex_alloc;
  while ((sequences + 1) * sizeof(seqinfo_t) > seqindex_alloc)
    {
      seqindex_alloc += memchunk;
    }
  if (seqindex_alloc > seqindex_alloc_old)
    {
      seqindex = (seqinfo_t *) util.xrealloc(seqindex, seqindex_alloc);
    }

  /* update index */
  seqinfo_t * seqindex_p = seqindex + sequences;
  seqindex_p->headerlen = headerlength;
  seqindex_p->seqlen = sequencelength;
  seqindex_p->header_p = header_p;
  seqindex_p->seq_p = sequence_p;
  seqindex_p->qual_p = quality_p;
  seqindex_p->size = abundance;

  /* update statistics */
  ++sequences;
  nucleotides += sequencelength;
  longest = std::max((uint64_t)sequencelength, longest);
  shortest = std::min((uint64_t)sequencelength, shortest);
  longestheader = std::max((uint64_t)headerlength, longestheader);
}
/******************************************************************************/
auto Vsearch_Database::read(std::vector<std::string>& names,
             std::vector<std::string>& seqs,
             std::vector<std::int64_t>& abundances) -> void
{

    longest = 0;
    shortest = std::numeric_limits<uint64_t>::max();  // refactoring: direct initialization
    longestheader = 0;
    sequences = 0;
    nucleotides = 0;

    /* allocate space for data */
    dataalloc = 0;
    datap = nullptr;
    datalen = 0;

    /* allocate space for index */
    seqindex_alloc = 0;
    seqindex = nullptr;

    for (int i = 0; i < seqs.size(); i++) {
        size_t const sequencelength = seqs[i].length();
        int64_t const abundance = abundances[i];

        add(names[i].c_str(), seqs[i].c_str(), nullptr,
            names[i].length(), sequencelength, abundance);
    }
}
/******************************************************************************/

