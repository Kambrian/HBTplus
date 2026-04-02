/*
 * Copyright (C) 2007-2026 Jiaxin Han and contributors
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef MYSORT_HEADER_INCLUDED
#define MYSORT_HEADER_INCLUDED

#ifdef _OPENMP
#include "parallel_stable_sort.h"
//the parallelization is done internally. so it's preferred to call the parallelsort outside an omp region.
#define MYSORT pss::parallel_stable_sort
#else
#define MYSORT sort
#endif

#endif