#ifndef VSEARCH_DBINDEX_h
#define VSEARCH_DBINDEX_h

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

#include "bitmap.h"
#include "db.h"
#include "unique.h"
#include <cstdio>  // std::FILE
#include <cstdint>  // uint64_t

constexpr unsigned int bitmap_threshold = 8;

struct uhandle_s;

/******************************************************************************/
class Vsearch_DBIndex {

public:

    Vsearch_DBIndex(Vsearch_Database* d,
                    Vsearch_Bitmap* b,
                    Vsearch_Unique* u) {
        db = d;
        bitmap = b;
        unique = u;
    }
    ~Vsearch_DBIndex() { free(); }

    auto prepare() -> void;

    auto addsequence(unsigned int seqno) -> void;

    auto addallsequences() -> void;

    auto getbitmap(unsigned int kmer) -> unsigned char *;

    auto getmatchcount(unsigned int kmer) -> unsigned int;

    auto getmatchlist(unsigned int kmer) -> unsigned int *;

    auto getmapping(unsigned int index) -> unsigned int;

    auto getcount() -> unsigned int;

    auto setWordLength(int64_t opt_wl) -> void { opt_wordlength = opt_wl; }

private:

    unsigned int * kmercount; /* number of matching seqnos for each kmer */
    uint64_t * kmerhash;  /* index into the list below for each kmer */
    unsigned int * kmerindex; /* the list of matching seqnos for kmers */
    struct bitmap_s * * kmerbitmap;
    unsigned int * dbindex_map;
    unsigned int dbindex_count;
    unsigned int kmerhashsize;
    uint64_t kmerindexsize;
    uhandle_s * dbindex_uh;
    unsigned int bitmap_mincount;
    int opt_wordlength;

    Vsearch_Database* db;
    Vsearch_Bitmap* bitmap;
    Vsearch_Unique* unique;
    Vsearch_Util util;

    void free() {
        util.xfree(kmerhash);
        util.xfree(kmerindex);
        util.xfree(kmercount);
        util.xfree(dbindex_map);

        for (auto kmer = 0U; kmer < kmerhashsize; kmer++)
        {
            if (kmerbitmap[kmer] != nullptr)
            {
                bitmap->free(kmerbitmap[kmer]);
            }
        }
        util.xfree(kmerbitmap);
        unique->exit(dbindex_uh);
    }
};
/******************************************************************************/

#endif
