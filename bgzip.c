/* bgzip.c -- Block compression/decompression utility.

   Copyright (C) 2008, 2009 Broad Institute / Massachusetts Institute of Technology
   Copyright (C) 2010, 2013-2019 Genome Research Ltd.

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notices and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#include <config.h>

#include <stdlib.h>
#include <string.h>
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <stdarg.h>
#include <inttypes.h>
#include "htslib/bgzf.h"
#include "htslib/hts.h"

#define isatty(file) FALSE

#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#endif

static const int WINDOW_SIZE = 64 * 1024;

static void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(EXIT_FAILURE);
}

static int ask_yn()
{
    char line[1024];
    if (fgets(line, sizeof line, stdin) == NULL)
        return 0;
    return line[0] == 'Y' || line[0] == 'y';
}

static int confirm_overwrite(const char *fn)
{
    int save_errno = errno;
    int ret = 0;

    if (isatty(STDIN_FILENO)) {
        fprintf(stderr, "[bgzip] %s already exists; do you wish to overwrite (y or n)? ", fn);
        if (ask_yn()) ret = 1;
    }

    errno = save_errno;
    return ret;
}

static int known_extension(const char *ext)
{
    static const char *known[] = {
        "gz", "bgz", "bgzf",
        NULL
    };

    const char **p;
    for (p = known; *p; p++)
        if (_stricmp(ext, *p) == 0) return 1;
    return 0;
}

static int confirm_filename(int *is_forced, const char *name, const char *ext)
{
    if (*is_forced) {
        (*is_forced)--;
        return 1;
    }

    if (!isatty(STDIN_FILENO))
        return 0;

    fprintf(stderr, "[bgzip] .%s is not a known extension; do you wish to decompress to %s (y or n)? ", ext, name);
    return ask_yn();
}