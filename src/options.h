#ifndef OPTIONS_H_
#define OPTIONS_H_

/*
 *  Modified by Sarah Westcott on 2/26/26.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 * This class will store the global options for the vsearch code
 *
 */

#include <limits>

constexpr int64_t default_fasta_width = 80;
constexpr int64_t default_fastq_tail = 4;
constexpr int64_t default_maxseqlength = 50000;
constexpr int64_t default_ascii_offset = 33;
constexpr char alternative_ascii_offset = 64;
constexpr auto dbl_max = std::numeric_limits<double>::max();
static constexpr auto int_max = std::numeric_limits<int>::max();
static constexpr auto long_min = std::numeric_limits<long>::min();
constexpr int64_t default_max_quality = 41;
constexpr auto int64_max = std::numeric_limits<int64_t>::max();
std::string const default_quality_padding = "IIIIIIII";  // Q40 with an offset of 33
std::string const alternative_quality_padding = "hhhhhhhh";  // Q40 with an offset of 64
std::string const default_sequence_padding = "NNNNNNNN";


class Vsearch_Options {

public:

    static Vsearch_Options* getInstance() {
        if( _uniqueInstance == 0) {
            _uniqueInstance = new Vsearch_Options();
        }
        return _uniqueInstance;
    }

    void reset() {
        opt_abskew = 2.0;
        opt_dn = 1.4;
        opt_xn = 8.0;
        opt_mindiv = 0.8;
        opt_minh = 0.28;
        opt_threads = 1;
        opt_chimeras_parts = 0; // will use 1 part per 100bp

        opt_maxsizeratio = dbl_max;
        opt_id = -1.0;
        opt_mindiffs = 3;
        opt_maxaccepts = 1;
        opt_maxrejects = -1;
        opt_self = 0;
        opt_selfid = 0;
        opt_strand = 1;
        opt_iddef = 2;
        opt_chimeras_parents_max = 3;
        opt_chimeras_diff_pct = 0.0;
        opt_chimeras_length_min = 10;
        //opt_n_mismatch = false;
        opt_mismatch = -4;
        opt_match = 2;
        opt_gap_extension_query_interior = 2;
        opt_gap_extension_query_left = 1;
        opt_gap_extension_query_right = 1;
        opt_gap_extension_target_interior = 2;
        opt_gap_extension_target_left = 1;
        opt_gap_extension_target_right = 1;
        opt_gap_open_query_interior = 20;
        opt_gap_open_query_left = 2;
        opt_gap_open_query_right = 2;
        opt_gap_open_target_interior = 20;
        opt_gap_open_target_left = 2;
        opt_gap_open_target_right = 2;
        opt_minwordmatches = -1;
        opt_maxqsize = int_max;
        opt_mintsize = 0;
        opt_minsizeratio = 0.0;
        opt_minqt = 0.0;
        opt_minsl = 0.0;
        opt_maxqt = dbl_max;
        opt_maxsl = dbl_max;
        opt_idprefix = 0;
        opt_idsuffix = 0;
        opt_leftjust = 0;
        opt_weak_id = -1.0;
        opt_maxsubs = int_max;
        opt_maxgaps = int_max;
        opt_maxdiffs = int_max;
        opt_rightjust = 0;
        opt_mincols = 0;
        opt_query_cov = 0.0;
        opt_target_cov = 0.0;
        opt_maxid = 1.0;
        opt_mid = 0.0;
        opt_unoise_alpha = 2.0;
        neon_present = 0;
        ssse3_present = 0;
        sse2_present = 0;
        altivec_present = 0;
        mmx_present = 0;
        sse_present = 0;
        sse3_present = 0;
        sse41_present = 0;
        sse42_present = 0;
        popcnt_present = 0;
        avx_present = 0;
        avx2_present = 0;
    }

    // user setable
    double opt_abskew, opt_dn, opt_xn, opt_mindiv, opt_minh;

    double opt_id;
    double opt_maxsizeratio;
    int opt_mindiffs;
    int64_t opt_maxaccepts;
    int64_t opt_maxrejects;
    int64_t opt_threads;
    int64_t opt_self;
    int64_t opt_selfid;
    int64_t opt_strand;
    //bool opt_n_mismatch;
    int64_t opt_mismatch;
    int64_t opt_match;
    int opt_chimeras_parts;
    int opt_chimeras_parents_max;
    double opt_chimeras_diff_pct;
    int opt_chimeras_length_min;
    int opt_gap_extension_query_interior;
    int opt_gap_extension_query_left;
    int opt_gap_extension_query_right;
    int opt_gap_extension_target_interior;
    int opt_gap_extension_target_left;
    int opt_gap_extension_target_right;
    int opt_gap_open_query_interior;
    int opt_gap_open_query_left;
    int opt_gap_open_query_right;
    int opt_gap_open_target_interior;
    int opt_gap_open_target_left;
    int opt_gap_open_target_right;
    int64_t opt_iddef;
    int64_t opt_minwordmatches;
    int64_t opt_maxqsize;
    int64_t opt_mintsize;
    double opt_minsizeratio;
    double opt_minqt;
    double opt_minsl;
    double opt_maxqt;
    double opt_maxsl;
    int64_t opt_idprefix;
    int64_t opt_idsuffix;
    int64_t opt_leftjust;
    double opt_weak_id;
    int64_t opt_maxsubs;
    int64_t opt_maxgaps;
    int64_t opt_rightjust;
    int64_t opt_mincols;
    double opt_query_cov;
    double opt_target_cov;
    double opt_maxid;
    double opt_mid;
    int64_t opt_maxdiffs;
    int64_t neon_present;
    double opt_unoise_alpha;
    int64_t ssse3_present;
    int64_t sse2_present;
    int64_t altivec_present;
    int64_t mmx_present;
    int64_t sse_present;
    int64_t sse3_present;
    int64_t sse41_present;
    int64_t sse42_present;
    int64_t popcnt_present;
    int64_t avx_present;
    int64_t avx2_present;

private:

    static Vsearch_Options* _uniqueInstance;
    Vsearch_Options( const Vsearch_Options& ); // Disable copy constructor
    void operator=( const Vsearch_Options& ); // Disable assignment operator
    Vsearch_Options() {
       reset();
   }
};

#endif
