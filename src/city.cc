// Copyright (c) 2011 Google, Inc.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
// CityHash, by Geoff Pike and Jyrki Alakuijala
//
// This file provides CityHash64() and related functions.
//
// It's probably possible to create even faster hash functions by
// writing a program that systematically explores some of the space of
// possible hash functions, by using SIMD instructions, or by
// compromising on hash quality.

//#include "config.h"

#include "city.h"
#include "os_byteswap.h"

#include <algorithm>  // std::swap
#include <cassert>
#include <cstdint>  // int32_t, uint8_t, uint32_t, uint64_t
#include <cstdio>  // std::size_t
#include <cstring>  // std::memcpy, std::memset
#include <utility> // std::pair, std::make_pair

// Some primes between 2^63 and 2^64 for various uses.
constexpr uint64_t k0 = 0xc3a5c85c97cb3127ULL;  // uint128 seeds?
constexpr uint64_t k1 = 0xb492b66fbe98f273ULL;  // uint128 seeds?
constexpr uint64_t k2 = 0x9ae16a3b2f90404fULL;  // uint64 initialization value?

// shift values
constexpr auto eighteen = 18U;
constexpr auto twentyone = 21U;
constexpr auto twentyfive = 25U;
//constexpr auto twentyseven = 27U;
constexpr auto thirty = 30U;
constexpr auto thirtythree = 33U;
//constexpr auto thirtyfive = 35U;
constexpr auto thirtyseven = 37U;
constexpr auto fortytwo = 42U;
constexpr auto fortythree = 43U;
constexpr auto fortyfour = 44U;
constexpr auto fortyseven = 47U;
//constexpr auto fortynine = 49U;
//constexpr auto fiftythree = 53U;
constexpr auto sixtyfour = 64U;

constexpr auto uint32_in_expected_order(uint32_t src) noexcept -> uint32_t {
#ifdef WORDS_BIGENDIAN
    return bswap_32(src);
#else
    return src;
#endif
}
constexpr auto uint64_in_expected_order(uint64_t src) noexcept -> uint64_t {
#ifdef WORDS_BIGENDIAN
    return bswap_64(src);
#else
    return src;
#endif
}
constexpr auto likely(bool condition) noexcept -> bool {
#if HAVE_BUILTIN_EXPECT
    static constexpr auto is_expected_to_be_true = 1L;
    // !!condition converts to a strict boolean value
    return __builtin_expect(!!condition, is_expected_to_be_true);
#else
    return condition;
#endif
}

/******************************************************************************/
auto Vsearch_CityHash::Uint128Low64(const uint128 & a_pair) -> uint64_t {
    return a_pair.first;
  }
/******************************************************************************/
auto Vsearch_CityHash::Uint128High64(const uint128 & a_pair) -> uint64_t {
    return a_pair.second;
  }
/******************************************************************************/
auto Vsearch_CityHash::unaligned_load64(const char * src) -> uint64_t {
    uint64_t result = 0;
    std::memcpy(&result, src, sizeof(result));
    return result;
  }
/******************************************************************************/
auto Vsearch_CityHash::unaligned_load32(const char * src) -> uint32_t {
    uint32_t result = 0;
    std::memcpy(&result, src, sizeof(result));
    return result;
  }
/******************************************************************************/
auto Vsearch_CityHash::Fetch64(const char * src) -> uint64_t {
    return uint64_in_expected_order(unaligned_load64(src));
  }
/******************************************************************************/
auto Vsearch_CityHash::Fetch32(const char * src) -> uint32_t {
    return uint32_in_expected_order(unaligned_load32(src));
  }
/******************************************************************************/
  // Bitwise right rotate.  Normally this will compile to a single
  // instruction, especially if the shift is a manifest constant.
  // C++20 refactoring: std::rotr()
auto Vsearch_CityHash::Rotate(uint64_t const val, unsigned int const shift) -> uint64_t {
    assert(shift != sixtyfour);  // shifting by 64 yields undefined results
    return shift == 0 ? val : ((val >> shift) | (val << (sixtyfour - shift)));
  }
/******************************************************************************/
auto Vsearch_CityHash::ShiftMix(uint64_t const val) -> uint64_t {
    return val ^ (val >> fortyseven);
  }
/******************************************************************************/
// Hash 128 input bits down to 64 bits of output.
// This is intended to be a reasonably good hash function.
auto Vsearch_CityHash::Hash128to64(const uint128 & a_pair) -> uint64_t {
    // Murmur-inspired hashing.
    static constexpr auto divider = 47U;
    static constexpr uint64_t kMul = 0x9ddfea08eb382d69ULL;
    auto reg_a = (Uint128Low64(a_pair) ^ Uint128High64(a_pair)) * kMul;
    reg_a ^= (reg_a >> divider);
    auto reg_b = (Uint128High64(a_pair) ^ reg_a) * kMul;
    reg_b ^= (reg_b >> divider);
    reg_b *= kMul;
    return reg_b;
  }
/******************************************************************************/
auto Vsearch_CityHash::HashLen16(uint64_t const first, uint64_t const second) -> uint64_t {
    return Hash128to64(uint128(first, second));
}
/******************************************************************************/
auto Vsearch_CityHash::HashLen16(uint64_t const first, uint64_t const second, uint64_t const mul) -> uint64_t {
    // Murmur-inspired hashing.
    auto reg_a = (first ^ second) * mul;
    reg_a ^= (reg_a >> fortyseven);
    auto reg_b = (second ^ reg_a) * mul;
    reg_b ^= (reg_b >> fortyseven);
    reg_b *= mul;
    return reg_b;
  }
/******************************************************************************/
auto Vsearch_CityHash::HashLen0to16(const char * seq, std::size_t const len) -> uint64_t {
    if (len >= 8) {
      const uint64_t mul = k2 + (len * 2);
      const uint64_t a = Fetch64(seq) + k2;
      const uint64_t b = Fetch64(seq + len - 8);
      const uint64_t c = (Rotate(b, thirtyseven) * mul) + a;
      const uint64_t d = (Rotate(a, twentyfive) + b) * mul;
      return HashLen16(c, d, mul);
    }
    if (len >= 4) {
      const uint64_t mul = k2 + (len * 2);
      const uint64_t a = Fetch32(seq);
      return HashLen16(len + (a << 3U), Fetch32(seq + len - 4), mul);
    }
    if (len > 0) {
      uint8_t const a = seq[0];
      uint8_t const b = seq[len >> 1U];
      uint8_t const c = seq[len - 1];
      const uint32_t y = static_cast<uint32_t>(a) + (static_cast<uint32_t>(b) << 8U);
      const uint32_t z = len + (static_cast<uint32_t>(c) << 2U);
      return ShiftMix((y * k2) ^ (z * k0)) * k2;
    }
    return k2;  // initialization value for empty sequences
  }
/******************************************************************************/
// This probably works well for 16-byte strings as well, but it may be overkill
// in that case.
auto Vsearch_CityHash::HashLen17to32(const char * seq, std::size_t const len) -> uint64_t {
    const uint64_t mul = k2 + (len * 2);
    const uint64_t a = Fetch64(seq) * k1;
    const uint64_t b = Fetch64(seq + 8);
    const uint64_t c = Fetch64(seq + len - 8) * mul;
    const uint64_t d = Fetch64(seq + len - 16) * k2;
    return HashLen16(Rotate(a + b, fortythree) + Rotate(c, thirty) + d,
                     a + Rotate(b + k2, eighteen) + c, mul);
  }
/******************************************************************************/
// Return a 16-byte hash for 48 bytes.  Quick and dirty.
// Callers do best to use "random-looking" values for a and b.
auto Vsearch_CityHash::WeakHashLen32WithSeeds(uint64_t const w, uint64_t const x,
                                              uint64_t const y, uint64_t const z,
                                              uint64_t a, uint64_t b) -> std::pair<uint64_t, uint64_t> {
    a += w;
    b = Rotate(b + a + z, twentyone);
    const uint64_t c = a;
    a += x;
    a += y;
    b += Rotate(a, fortyfour);
    return std::make_pair(a + z, b + c);
  }
/******************************************************************************/
// Return a 16-byte hash for s[0] ... s[31], a, and b.  Quick and dirty.
auto Vsearch_CityHash::WeakHashLen32WithSeeds(const char * seq, uint64_t a, uint64_t b)
    -> std::pair<uint64_t, uint64_t> {
    return WeakHashLen32WithSeeds(Fetch64(seq),
                                  Fetch64(seq + 8),
                                  Fetch64(seq + 16),
                                  Fetch64(seq + 24),
                                  a,
                                  b);
  }
/******************************************************************************/
// Return an 8-byte hash for 33 to 64 bytes.
auto Vsearch_CityHash::HashLen33to64(const char * seq, std::size_t len) -> uint64_t {
    const uint64_t mul = k2 + (len * 2);
    uint64_t a = Fetch64(seq) * k2;
    uint64_t b = Fetch64(seq + 8);
    const uint64_t c = Fetch64(seq + len - 24);
    const uint64_t d = Fetch64(seq + len - 32);
    const uint64_t e = Fetch64(seq + 16) * k2;
    const uint64_t f = Fetch64(seq + 24) * 9;
    const uint64_t g = Fetch64(seq + len - 8);
    const uint64_t h = Fetch64(seq + len - 16) * mul;
    const uint64_t u = Rotate(a + g, fortythree) + ((Rotate(b, thirty) + c) * 9);
    const uint64_t v = ((a + g) ^ d) + f + 1;
    const uint64_t w = byteSwap.bswap_64((u + v) * mul) + h;
    const uint64_t x = Rotate(e + f, fortytwo) + c;
    const uint64_t y = (byteSwap.bswap_64((v + w) * mul) + g) * mul;
    const uint64_t z = e + f + c;
    a = byteSwap.bswap_64(((x + z) * mul) + y) + b;
    b = ShiftMix(((z + a) * mul) + d + h) * mul;
    return b + x;
  }
/******************************************************************************/
auto Vsearch_CityHash::CityHash64(const char * seq, std::size_t len) -> uint64_t {
  if (len <= 16) {
    return HashLen0to16(seq, len);
  }
  if (len <= 32) {
    return HashLen17to32(seq, len);
  }
  if (len <= 64) {
    return HashLen33to64(seq, len);
  }

  // For strings over 64 bytes we hash the end first, and then as we
  // loop we keep 56 bytes of state: v, w, x, y, and z.
  uint64_t x = Fetch64(seq + len - 40);
  uint64_t y = Fetch64(seq + len - 16) + Fetch64(seq + len - 56);
  uint64_t z = HashLen16(Fetch64(seq + len - 48) + len, Fetch64(seq + len - 24));
  auto v = WeakHashLen32WithSeeds(seq + len - 64, len, z);
  auto w = WeakHashLen32WithSeeds(seq + len - 32, y + k1, x);
  x = (x * k1) + Fetch64(seq);

  // Decrease len to the nearest multiple of 64, and operate on 64-byte chunks.
  len = (len - 1) & ~static_cast<std::size_t>(63);
  do {
    x = Rotate(x + y + v.first + Fetch64(seq + 8), thirtyseven) * k1;
    y = Rotate(y + v.second + Fetch64(seq + 48), fortytwo) * k1;
    x ^= w.second;
    y += v.first + Fetch64(seq + 40);
    z = Rotate(z + w.first, thirtythree) * k1;
    v = WeakHashLen32WithSeeds(seq, v.second * k1, x + w.first);
    w = WeakHashLen32WithSeeds(seq + 32, z + w.second, y + Fetch64(seq + 16));
    std::swap(z, x);
    seq += 64;
    len -= 64;
  } while (len != 0);
  return HashLen16(HashLen16(v.first, w.first) + (ShiftMix(y) * k1) + z,
                   HashLen16(v.second, w.second) + x);
}
/******************************************************************************/
