#ifndef VSEARCH_DATABASE_h
#define VSEARCH_DATABASE_h

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

#include <cstdint>  // uint64_t
#include <cstdio>  // std::size_t
#include <string>
#include <vector>

#include "options.h"
#include "util.h"

/******************************************************************************/
struct seqinfo_s
{
  std::size_t header_p;
  std::size_t seq_p;
  std::size_t qual_p;
  unsigned int headerlen;
  unsigned int seqlen;
  unsigned int size;
};
/******************************************************************************/
using seqinfo_t = struct seqinfo_s;
/******************************************************************************/
class Vsearch_Database {

public:

    Vsearch_Database() {
        sequences = 0;
        nucleotides = 0;
        longest = 0;
        shortest = 0;
        longestheader = 0;

        dataalloc = 0;
        datalen = 0;
        seqindex_alloc = 0;

        seqindex = nullptr;
        datap = nullptr;
        opts = Vsearch_Options::getInstance();
    }
    ~Vsearch_Database() { free(); }

    auto getheader(uint64_t seqno) -> char * {
        return datap + seqindex[seqno].header_p; }
    auto getsequence(uint64_t seqno) -> char * {
        return datap + seqindex[seqno].seq_p; }
    auto getabundance(uint64_t seqno) -> uint64_t {
        return seqindex[seqno].size; }
    auto getsequencelen(uint64_t seqno) -> uint64_t {
        return seqindex[seqno].seqlen; }
    auto getheaderlen(uint64_t seqno) -> uint64_t {
        return seqindex[seqno].headerlen; }
    auto getsequencecount() -> uint64_t {
        return sequences; }
    auto getnucleotidecount() -> uint64_t {
        return nucleotides; }
    auto getlongestheader() -> uint64_t {
        return longestheader; }
    auto getlongestsequence() -> uint64_t {
        return longest; }
    auto getshortestsequence() -> uint64_t {
        return shortest; }

    auto read(std::vector<std::string>& names,
                 std::vector<std::string>& seqs,
                 std::vector<std::int64_t>& abundances) -> void;

private:

    Vsearch_Util util;
    Vsearch_Options* opts;
    char * datap;
    seqinfo_t * seqindex;

    uint64_t sequences, nucleotides, longest, shortest, longestheader;
    uint64_t dataalloc, datalen;
    size_t seqindex_alloc;

    auto add(char const * header,
             char const * sequence,
             char const * quality,
             size_t const headerlength,
             size_t const sequencelength,
             int64_t const abundance) -> void;

    void free() {
        if (datap != nullptr)
        {
            util.xfree(datap);
        }
        if (seqindex != nullptr)
        {
            util.xfree(seqindex);
        }
    }
};
/******************************************************************************/

#endif

