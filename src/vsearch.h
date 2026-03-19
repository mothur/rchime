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

#define _GNU_SOURCE 1
// #define __STDC_CONSTANT_MACROS 1
// #define __STDC_FORMAT_MACROS 1
// #define __STDC_LIMIT_MACROS 1
#define __restrict

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "Rcpp.h"
#include <cctype>
#include <cfloat>
#include <cinttypes>
#include <clocale>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#include <cstdio>  // std::size_t
#include <cstdlib>
#include <cstring>
#include <ctime>  // replace with std::chrono
#include <map>
#include <set>
#include <string>
#include <cassert>
#include <limits>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>  // refactoring: redundant with <ctime>?
#include <pthread.h>
#include <getopt.h>
#include <fcntl.h>
#include <unistd.h>

#define PROG_NAME PACKAGE
#define PROG_VERSION PACKAGE_VERSION

#ifdef __x86_64__

#define PROG_CPU "x86_64"
#include <x86intrin.h>

#elif __PPC__

#ifdef __LITTLE_ENDIAN__
#define PROG_CPU "ppc64le"
#include <altivec.h>
#undef bool
#else
#error Big endian ppc64 CPUs not supported
#endif

#elif __aarch64__

#define PROG_CPU "aarch64"
#include <arm_neon.h>

#else

#define PROG_CPU "simde"
#define SIMDE_ENABLE_NATIVE_ALIASES
#include <simde/x86/avx512.h>
#endif


#ifdef _WIN32

#define PROG_OS "win"
#include <windows.h>
#include <psapi.h>
#include <shlwapi.h>

#elif __APPLE__

#define PROG_OS "macos"
#include <sys/sysctl.h>
#include <sys/resource.h>

#elif __linux__

#define PROG_OS "linux"
#include <sys/sysinfo.h>
#include <sys/resource.h>

#elif __FreeBSD__

#define PROG_OS "freebsd"
#include <sys/sysinfo.h>
#include <sys/resource.h>

#elif __NetBSD__

#define PROG_OS "netbsd"
#include <sys/resource.h>
#include <sys/types.h>
/* Alters behavior, but NetBSD 7 does not have getopt_long_only() */
#define getopt_long_only getopt_long

#else

#define PROG_OS "unknown"
#include <sys/sysinfo.h>
#include <sys/resource.h>

#endif


#define PROG_ARCH PROG_OS "_" PROG_CPU

#ifdef HAVE_DLFCN_H
#include <dlfcn.h>
#endif

#ifdef HAVE_ZLIB_H
#include <zlib.h>
#endif

#ifdef HAVE_BZLIB_H
#include <bzlib.h>
#endif

