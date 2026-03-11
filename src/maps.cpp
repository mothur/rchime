/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2025, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
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

#include "maps.h"

/******************************************************************************/
auto Vsearch_Maps::map_uppercase(char const nucleotide) -> char {
  return static_cast<char>(chrmap_upcase_vector[to_uchar(nucleotide)]);
}

/******************************************************************************/
auto Vsearch_Maps::map_2bit(char const nucleotide) -> unsigned int {
  return chrmap_2bit[to_uchar(nucleotide)];
}

/******************************************************************************/
auto Vsearch_Maps::map_4bit(char const nucleotide) -> unsigned char {
  return chrmap_4bit[to_uchar(nucleotide)];
}

/******************************************************************************/
auto Vsearch_Maps::map_complement(char const nucleotide) -> char {
  return static_cast<char>(chrmap_complement[to_uchar(nucleotide)]);
}

/******************************************************************************/
auto Vsearch_Maps::map_mask_ambig(char const nucleotide) -> unsigned int {
  return chrmap_mask_ambig[to_uchar(nucleotide)];
}

/******************************************************************************/
auto Vsearch_Maps::map_mask_lower(char const nucleotide) -> unsigned int {
  return chrmap_mask_lower[to_uchar(nucleotide)];
}

/******************************************************************************/
auto Vsearch_Maps::is_equivalent_4bit(char const lhs, char const rhs) -> bool {
  auto const lhs_unsigned = map_4bit(lhs);
  auto const rhs_unsigned = map_4bit(rhs);
  return ((lhs_unsigned & rhs_unsigned) != 0);
}

/******************************************************************************/
auto Vsearch_Maps::is_equivalent_4bit_rhs(char const lhs, char const rhs) -> bool {
  auto const lhs_unsigned = to_uchar(lhs);
  auto const rhs_unsigned = map_4bit(rhs);
  return ((lhs_unsigned & rhs_unsigned) != 0);
}

/******************************************************************************/
auto Vsearch_Maps::is_ambiguous_4bit(unsigned char const nucleotide) -> bool {
  return chrmap_ambiguous_4bit[nucleotide];
}

/******************************************************************************/
auto Vsearch_Maps::is_same_4bit(char const lhs, char const rhs) -> bool {
  return map_4bit(lhs) == map_4bit(rhs);
}
/******************************************************************************/
