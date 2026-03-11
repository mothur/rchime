#ifndef VSEARCH_ALIGNSIMD_h
#define VSEARCH_ALIGNSIMD_h

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

#include "db.h"
#include "util.h"
#include "maps.h"
#include "vsearch.h"

using CELL = signed short;
using WORD = unsigned short;
using BYTE = unsigned char;

struct s16info_s;

constexpr auto matrix_size = 16;
constexpr auto CHANNELS = 8;
constexpr auto CDEPTH = 4;
constexpr auto maxseqlenproduct = 25000000LL;

#ifdef __PPC__

using VECTOR_SHORT = __vector signed short;

#elif defined __aarch64__

using VECTOR_SHORT = int16x8_t;

#elif defined(__x86_64__) || defined(SIMDE_VERSION)

using VECTOR_SHORT = __m128i;

#else

#error Unknown Architecture

#endif

/******************************************************************************/
class Vsearch_AlignSimd {

public:

    Vsearch_AlignSimd() = default;
    ~Vsearch_AlignSimd() = default;

    auto search16_init(CELL score_match,
                       CELL score_mismatch,
                       CELL penalty_gap_open_query_left,
                       CELL penalty_gap_open_target_left,
                       CELL penalty_gap_open_query_interior,
                       CELL penalty_gap_open_target_interior,
                       CELL penalty_gap_open_query_right,
                       CELL penalty_gap_open_target_right,
                       CELL penalty_gap_extension_query_left,
                       CELL penalty_gap_extension_target_left,
                       CELL penalty_gap_extension_query_interior,
                       CELL penalty_gap_extension_target_interior,
                       CELL penalty_gap_extension_query_right,
                       CELL penalty_gap_extension_target_right,
                       Vsearch_Options* opts) -> struct s16info_s *;


    auto search16_exit(s16info_s * searchinfo) -> void;


    auto search16_qprep(s16info_s * searchinfo, char * qseq, int qlen) -> void;


    auto search16(s16info_s * searchinfo,
                  unsigned int sequences,
                  unsigned int * seqnos,
                  CELL * pscores,
                  unsigned short * paligned,
                  unsigned short * pmatches,
                  unsigned short * pmismatches,
                  unsigned short * pgaps,
                  char * * pcigar,
                  Vsearch_Database* db, bool opt_n_mismatch) -> void;

private:

    Vsearch_Util util;
    Vsearch_Maps maps;

    std::array<std::array<int64_t, matrix_size>, matrix_size> scorematrix {{}};


    auto dprofile_fill16(CELL * dprofile_word,
                         CELL * score_matrix_word,
                         BYTE const * dseq) -> void;
    auto compute_score_min(struct s16info_s const & alignment) -> short;

    auto backtrack16(s16info_s * s,
                     char * dseq,
                     uint64_t dlen,
                     uint64_t offset,
                     uint64_t channel,
                     unsigned short * paligned,
                     unsigned short * pmatches,
                     unsigned short * pmismatches,
                     unsigned short * pgaps,
                     bool) -> void;
    auto pushop(s16info_s * s, char const newop) -> void;
    auto finishop(s16info_s * s) -> void;
    auto aligncolumns_rest(VECTOR_SHORT * Sm,
                           VECTOR_SHORT * hep,
                                              VECTOR_SHORT ** qp,
                                              VECTOR_SHORT QR_q_i,
                                              VECTOR_SHORT R_q_i,
                                              VECTOR_SHORT QR_q_r,
                                              VECTOR_SHORT R_q_r,
                                              VECTOR_SHORT QR_t_0,
                                              VECTOR_SHORT R_t_0,
                                              VECTOR_SHORT QR_t_1,
                                              VECTOR_SHORT R_t_1,
                                              VECTOR_SHORT QR_t_2,
                                              VECTOR_SHORT R_t_2,
                                              VECTOR_SHORT QR_t_3,
                                              VECTOR_SHORT R_t_3,
                                              VECTOR_SHORT h0,
                                              VECTOR_SHORT h1,
                                              VECTOR_SHORT h2,
                                              VECTOR_SHORT h3,
                                              VECTOR_SHORT f0,
                                              VECTOR_SHORT f1,
                                              VECTOR_SHORT f2,
                                              VECTOR_SHORT f3,
                                              VECTOR_SHORT * _h_min,
                                              VECTOR_SHORT * _h_max,
                                              int64_t ql,
                                              unsigned short * dir) -> void;
    auto aligncolumns_first(VECTOR_SHORT * Sm,
                            VECTOR_SHORT * hep,
                            VECTOR_SHORT ** qp,
                            VECTOR_SHORT QR_q_i,
                            VECTOR_SHORT R_q_i,
                            VECTOR_SHORT QR_q_r,
                            VECTOR_SHORT R_q_r,
                            VECTOR_SHORT QR_t_0,
                            VECTOR_SHORT R_t_0,
                            VECTOR_SHORT QR_t_1,
                            VECTOR_SHORT R_t_1,
                            VECTOR_SHORT QR_t_2,
                            VECTOR_SHORT R_t_2,
                            VECTOR_SHORT QR_t_3,
                            VECTOR_SHORT R_t_3,
                            VECTOR_SHORT h0,
                            VECTOR_SHORT h1,
                            VECTOR_SHORT h2,
                            VECTOR_SHORT h3,
                            VECTOR_SHORT f0,
                            VECTOR_SHORT f1,
                            VECTOR_SHORT f2,
                            VECTOR_SHORT f3,
                            VECTOR_SHORT * _h_min,
                            VECTOR_SHORT * _h_max,
                            VECTOR_SHORT Mm,
                            VECTOR_SHORT M_QR_t_left,
                            VECTOR_SHORT M_R_t_left,
                            VECTOR_SHORT M_QR_q_interior,
                            VECTOR_SHORT M_QR_q_right,
                            int64_t ql,
                            unsigned short * dir) -> void;
};
/******************************************************************************/

#endif
