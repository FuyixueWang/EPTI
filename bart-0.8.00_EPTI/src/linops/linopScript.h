/* Copyright 2014-2015. The Regents of the University of California.
 * Copyright 2016-2018. Martin Uecker.
 * All rights reserved. Use of this source code is governed by
 * a BSD-style license which can be found in the LICENSE file.
 */

#include <complex.h>

#ifndef __LINOPSCRIPT_H
#define __LINOPSCRIPT_H

#include "misc/cppwrap.h"
#include "misc/types.h"

#include "linops/linop.h"


#define MAX_LINOPS 20
#define MAX_LINSCRIPTS 20
#define MAX_FILES 10
#define MAX_TENSORS_BYSLICE 10

extern TYPEID linop_data_s;

typedef void (*lop_fun_t)(const linop_data_t* _data, complex float* dst, const complex float* src);
typedef void (*lop_p_fun_t)(const linop_data_t* _data, float lambda, complex float* dst, const complex float* src);
typedef void (*del_fun_t)(const linop_data_t* _data);

struct operator_s;
struct operator_p_s;

struct linop_s;

long getLinopOutCounter();
const struct linop_s** getLinopsOutVec();
void FreeOutLinops();
const struct linop_s* getLinopOut(long n);

long * getFdims(long i);
complex float* getDataFile(long i);

const struct linop_s* getLinopScriptFromFile(const char *FN, const complex float* InDimsC, long nInSets);
void FreeLinops(bool NotLast);
void ReadScriptFiles(char* argv[],long n);
void ReadScriptFiles_gpu(char* argv[],long n); // add
void ClearReadScriptFiles( char* argv[],long n);

#include "misc/cppwrap.h"

#endif // __LINOP_H

