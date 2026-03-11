#ifndef VSEARCH_MAIN_h
#define VSEARCH_MAIN_h

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
#include "options.h"
#include "db.h"
#include "dbindex.h"
#include "searchcore.h"
#include "cigar.h"
#include "maps.h"

constexpr auto maxparents = 20; /* max, could be fewer */
constexpr auto maxparts = 100;
constexpr auto window = 32;
constexpr auto few = 4;
constexpr auto maxcandidates = few * maxparts;
constexpr auto rejects = 16;
constexpr auto chimera_id = 0.55;
constexpr int dust_window = 64;

using CELL = signed short;
using WORD = unsigned short;
using BYTE = unsigned char;

/******************************************************************************/
/* information for each query sequence to be checked */
struct chimera_info_s
{
    int query_alloc = 0; /* the longest query sequence allocated memory for */
    int head_alloc = 0; /* the longest header allocated memory for */
    int part_alloc = 0; /* the longest query part allocated memory for */

    int query_no = 0;
    std::vector<char> query_head;
    int query_head_len = 0;
    int query_size = 0;
    std::vector<char> query_seq;
    int query_len = 0;

    std::array<struct searchinfo_s, maxparts> si {{}};

    std::array<unsigned int, maxcandidates> cand_list {{}};
    int cand_count = 0;

    struct s16info_s * s = nullptr;
    std::array<CELL, maxcandidates> snwscore {{}};
    std::array<unsigned short, maxcandidates> snwalignmentlength {{}};
    std::array<unsigned short, maxcandidates> snwmatches {{}};
    std::array<unsigned short, maxcandidates> snwmismatches {{}};
    std::array<unsigned short, maxcandidates> snwgaps {{}};
    std::array<int64_t, maxcandidates> nwscore {{}};
    std::array<int64_t, maxcandidates> nwalignmentlength {{}};
    std::array<int64_t, maxcandidates> nwmatches {{}};
    std::array<int64_t, maxcandidates> nwmismatches {{}};
    std::array<int64_t, maxcandidates> nwgaps {{}};
    std::vector<char *> nwcigar = std::vector<char *>(maxcandidates);  // this is a test

    int match_size = 0;
    std::vector<int> match;
    std::vector<int> insert;
    std::vector<int> smooth;
    std::vector<int> maxsmooth;

    std::vector<double> scan_p;
    std::vector<double> scan_q;

    int parents_found = 0;
    std::array<int, maxparents> best_parents {{}};
    std::array<int, maxparents> best_start {{}};
    std::array<int, maxparents> best_len {{}};

    int best_target = 0;
    char * best_cigar = nullptr;

    std::vector<int> maxi;  // longest insertion per position
    std::vector<std::vector<char>> paln;
    std::vector<char> qaln;
    std::vector<char> diffs;
    std::vector<char> votes;
    std::vector<char> model;
    std::vector<bool> ignore;

    double best_h = 0;
};
/******************************************************************************/
enum struct Status : unsigned char {
    no_parents,   // (0) non-chimeric
    no_alignment, // (1) score < 0, non-chimeric
    low_score,    // (2) score < minh, non-chimeric
    suspicious,   // (3) score >= minh, not available with uchime2_denovo and uchime3_denovo
    chimeric      // (4) score >= minh && divdiff >= opt_mindiv && ...
};
/******************************************************************************/
struct parents_info_s
{
    int cand = -1;
    int start = -1;
    int len = 0;
};

/******************************************************************************/
class Vsearch_Main {

public:

    Vsearch_Main();
    ~Vsearch_Main();

    auto vmain(std::vector<std::string>& sequenceNames,
               std::vector<std::string>& sequences,
               std::vector<std::string>& refNames,
               std::vector<std::string>& refSeqs,
               std::vector<float>& abunds,
               unsigned int start, unsigned int stop) -> std::vector<ChimeHit2>;

    std::set<std::string> chimeraNames;

private:

    bool uchimeDeNovo;

    std::vector<ChimeHit2> results;
    std::vector<std::string>* queryNames;
    std::vector<std::string>* querySeqs;
    std::vector<std::int64_t>* queryAbunds;

    int parts, tophits;
    unsigned int seqno, referenceThreadStop;

    int chimera_count, nonchimera_count, borderline_count, total_count;
    int64_t chimera_abundance, nonchimera_abundance, borderline_abundance, total_abundance = 0;

    chimera_info_s * cia;
    Vsearch_Options* opts;
    Vsearch_Util util;
    Vsearch_Database* db;
    Vsearch_DBIndex* dbindex;
    Vsearch_Unique* unique;
    Vsearch_Minheap* minheap;
    Vsearch_AlignSimd* align;
    Vsearch_Searchcore* search;
    Vsearch_Bitmap* bitmap;
    Vsearch_Cigar* cigar;
    Vsearch_Maps* maps;

    auto cpu_features_detect() -> void;
    auto realloc_arrays(struct chimera_info_s * chimera_info) -> void;
    auto reset_matches(struct chimera_info_s * a_chimera_info) -> void;
    auto find_matches(struct chimera_info_s * chimera_info) -> void;
    auto compare_positions(const void * a, const void * b) -> int;
    auto find_best_parents(struct chimera_info_s * ci) -> int;
    auto find_total_alignment_length(struct chimera_info_s const * chimera_info) -> int;
    auto fill_max_alignment_length(struct chimera_info_s * chimera_info) -> void;
    auto fill_alignment_parents(struct chimera_info_s * ci) -> void;
    auto fill_in_alignment_string_for_query(struct chimera_info_s * chimera_info) -> void;
    auto fill_in_model_string_for_query(struct chimera_info_s * chimera_info) -> void;
    auto count_matches_with_parents(struct chimera_info_s const * chimera_info,
                                    int const alignment_length) -> std::array<int, maxparents>;
    auto compute_global_similarities_with_parents(
            std::array<int, maxparents> const & match_counts,
            int const alignment_length) -> std::array<double, maxparents>;
    auto compute_diffs(struct chimera_info_s const * ci,
                       std::vector<unsigned char> const & psym,
                       unsigned char const qsym) -> char;
    auto eval_parents(struct chimera_info_s * ci) -> Status;
    auto query_init(struct searchinfo_s * search_info) -> void;
    auto query_exit(struct searchinfo_s * search_info) -> void;
    auto partition_query(struct chimera_info_s * chimera_info) -> void;
    auto chimera_init(struct chimera_info_s * ci) -> void;
    auto chimera_exit(struct chimera_info_s * ci) -> void;
    // main chimera checking done here
    auto chimera_core(struct chimera_info_s * ci) -> uint64_t;

    // from mask.h
    auto dust_all(uint64_t numSeqs) -> void;
    auto dust(char * seq, int len) -> void;
    auto wo(int len, const char *s, int *beg, int *end) -> int;

};
/******************************************************************************/

#endif // VSEARCH_MAIN_h
