/* Copyright 2013, 2016. The Regents of the University of California.
 * Copyright 2015. Martin Uecker.
 * All rights reserved. Use of this source code is governed by
 * a BSD-style license which can be found in the LICENSE file.
 *
 * Authors: 
 * 2012, 2015 Martin Uecker <martin.uecker@med.uni-goettingen.de>
 * 2016 Jon Tamir <jtamir@eecs.berkeley.edu>
 */

#include <stdbool.h>
#include <complex.h>

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "num/multind.h"
#include "num/flpmath.h"
#include "num/init.h"
#include "num/ops.h"
#include "num/iovec.h"

#include "num/casorati.h"


#include "misc/mmio.h"
#include "misc/misc.h"
#include "misc/opts.h"
#include "misc/debug.h"

#include "linops/linop.h"
#include "linops/linopScript.h"
#include "linops/fmac.h"
#include "linops/sampling.h"
#include "linops/someops.h"
#include "linops/sampling.h"
#include "linops/finite_diff.h"
#include "num/fft.h"
#include "misc/utils.h"

#include "noncart/nufft.h"
#include "noncart/nudft.h"
#include "noncart/precond.h"

#include "lowrank/lrthresh.h"

#include "iter/italgos.h"
#include "iter/iter.h"
#include "iter/prox.h"
#include "iter/admm.h"
#include "iter/vec.h"
#include "iter/niht.h"

#include "iter/iter2.h"

#include "iter/prox.h"
#include "iter/thresh.h"

#include "linops/grad.h"
#include "linops/sum.h"
#include "linops/waveop.h"

#include "wavelet/wavthresh.h"


#include "misc/mri.h"
#include "misc/utils.h"

#include "grecon/optreg.h"

#ifndef DIMS
#define DIMS 16
#endif

// Dictionary?
#define MAX_DICT_ITEMS 100
#define MAX_DICT_LEN 20
#define DICT_ERR -777
static char Dict[MAX_DICT_ITEMS][MAX_DICT_LEN];
static long DictVals[MAX_DICT_ITEMS];
static long DictSize=0;

unsigned char tolowerx(unsigned char ch) {
    if (ch >= 'A' && ch <= 'Z')
        ch = 'a' + (ch - 'A');
    return ch;
 }

int strcasecmp(const char *s1, const char *s2) {
    const unsigned char *us1 = (const u_char *)s1,
                        *us2 = (const u_char *)s2;

    while (tolowerx(*us1) == tolowerx(*us2++))
        if (*us1++ == '\0')
            return (0);
    return (tolowerx(*us1) - tolowerx(*--us2));
}

void AddToDict(const char* item, const long val) {
    strcpy((char*)Dict[DictSize],item);
    DictVals[DictSize++]=val;

    debug_printf(DP_INFO,"Added \"%s\" to dictionary with value %d\n", item,val);
}
long LookInDict(const char* item) {
    long i;
    for(i=0;i<DictSize;i++) {
        // debug_printf(DP_INFO,"Comparing %s and %s \n",item,Dict[i]);
        // if(strcasecmp(item,Dict[i])==0) {
        if(strcmp(item,Dict[i])==0) {
            // debug_printf(DP_INFO,"Match!\n");
            return DictVals[i];
        }
    }
    return DICT_ERR;
}
long LookInDictOrStr(const char* item) {
    debug_printf(DP_INFO,"LookInDictOrStr for %s\n", item);
  char *temp;
// //   errno = 0;
//   long value = strtol(item,&temp,10); // using base 10
//   if(temp == item) {debug_printf(DP_INFO,"temp==item);\n"); }
//   if(*temp != '\0') {debug_printf(DP_INFO,"*temp != \\0);\n"); }
//   debug_printf(DP_INFO,"value %d\n",value);
    long value=LookInDict(item);
    if(value==DICT_ERR) {
        value = strtol(item,&temp,10);
        debug_printf(DP_INFO,"LookInDictOrStr: Parsed to %d\n", value);
    } else {
        debug_printf(DP_INFO,"LookInDictOrStr: Found in dictionary and returned %d\n", value);
    }
//   if (temp == item || *temp != '\0') {
//         // ((*val == LONG_MIN || *val == LONG_MAX) && errno == ERANGE))
//       value=LookInDict(item);
//       debug_printf(DP_INFO,"LookInDictOrStr: Found in dictionary and returned %d\n", value);
//   } else {
//       debug_printf(DP_INFO,"LookInDictOrStr: Parsed to %d\n", value);
//   }
  return value;
}
// End dictionary


#define NEW_OP_POINTER  LinopsVec[LinopCounter]

#define ADD_OP Sop = linop_chain(Sop,LinopsVec[LinopCounter]);md_copy_dims(DIMS, CurDims, linop_codomain(LinopsVec[LinopCounter])->dims);debug_printf(DP_INFO,"OutDims: ");debug_print_dims(DP_INFO,DIMS,CurDims);LinopCounter++;

#define GET_NEXT_TOKEN token = strtok(NULL, " ,.\n");

complex float* dataFiles[MAX_FILES];
long Fdims[MAX_FILES][DIMS];

const struct linop_s* LinopsVec[MAX_LINOPS];
long LinopCounter=0;

const struct linop_s* LinopsOutVec[MAX_LINSCRIPTS];
long LinopOutCounter=0;

long DefOp=0;

long getLinopOutCounter() { return LinopOutCounter+1-DefOp; }
const struct linop_s** getLinopsOutVec() { return &LinopsOutVec[DefOp]; }
const struct linop_s* getLinopOut(long n) { return LinopsOutVec[DefOp+n]; }
const struct linop_s* getLinopOutX(long n) { return LinopsOutVec[n]; }


long * getFdims(long i) { return Fdims[i]; }
complex float* getDataFile(long i) { return dataFiles[i]; }

void FreeOutLinops() {
    // debug_printf(DP_INFO,"Freeing %ld Out linops\n",LinopOutCounter);
    for(long i=0;i<LinopOutCounter;i++) {
        if(LinopsOutVec[i]) {
            linop_free(LinopsOutVec[i]); }
    }
    // debug_printf(DP_INFO,"FreeOutLinops done\n");
}

void FreeLinops(bool NotLast) {
    for(long i=0;i<LinopCounter-1;i++) {
        if(LinopsVec[i]) {
            linop_free(LinopsVec[i]); }
    }
    if(!NotLast) {
        linop_free(LinopsVec[LinopCounter-1]);
    }
   // debug_printf(DP_INFO,"FreeLinops done\n");
}
void ReadScriptFiles_gpu(char* argv[],long n) {
    long i;
    debug_printf(DP_INFO,"Reading files\n");
    for(i=0;i<n;i++) {        
        complex float* tmp = load_cfl(argv[i], DIMS, Fdims[i]);
#ifdef USE_CUDA
        dataFiles[i]=md_gpu_move(DIMS, Fdims[i], tmp, CFL_SIZE);
#else
        dataFiles[i]=tmp;
#endif
        debug_printf(DP_INFO,"Reading %s: %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld\n",argv[i],Fdims[i][0],Fdims[i][1],Fdims[i][2],Fdims[i][3],Fdims[i][4],Fdims[i][5],Fdims[i][6],Fdims[i][7],Fdims[i][8],Fdims[i][9],Fdims[i][10],Fdims[i][11],Fdims[i][12],Fdims[i][13],Fdims[i][14],Fdims[i][15]);
    }
    debug_printf(DP_INFO,"Finished reading files\n");
}

void ReadScriptFiles(char* argv[],long n) {
    long i;
    debug_printf(DP_INFO,"Reading files\n");
    for(i=0;i<n;i++) {        
        dataFiles[i] = load_cfl(argv[i], DIMS, Fdims[i]);
    }
    debug_printf(DP_INFO,"Finished reading files\n");
}

void ClearReadScriptFiles( char* argv[],long n) {
    long i;
    debug_printf(DP_INFO,"Clearing files' memory\n");
    for(i=0;i<n;i++) {
        // debug_printf(DP_INFO,"Clearing %s\n",argv[i]);
        unmap_cfl(DIMS, Fdims[i], dataFiles[i]);
    }
    debug_printf(DP_INFO,"Finished Clearing files' memory\n");
}

const struct linop_s* getLinopScriptFromFile(const char *FN, const complex float* InDimsC, long nInSets) {
    LinopOutCounter=0;
    LinopCounter=0;
  //  const struct linop_s* Sop = linop_identity_create(DIMS, CurDims);

    const struct linop_s* SopBackup;
    
    FILE * fp;
	char * line = NULL;
    char *token;
	size_t len = 0;
	ssize_t read;
    long i,j;
    debug_printf(DP_INFO,"getLinopScriptFromFile start, %d sets\n",nInSets);
    long StartDims[MAX_LINSCRIPTS][DIMS];
    long *StartDimsp=StartDims;
    for(long i=0;i<nInSets;i++) {
        for(long j=0;j<DIMS;j++) {
        //    debug_printf(DP_WARN, "!!!!!!!!!!! sets %d  \n",InDimsC[i*DIMS+j]);
            StartDims[i][j]=crealf(InDimsC[i*DIMS+j]); } }


    for(long d=0;d<nInSets;d++) {
        debug_printf(DP_INFO,"StartDims[..][%d]: ",d);
        debug_print_dims(DP_INFO,DIMS,StartDims[d]);
    }
    
    long CurDims[DIMS];
    md_copy_dims(DIMS, CurDims, StartDims[0]);

    const struct linop_s* Sop = linop_identity_create(DIMS, CurDims);

    fp = fopen(FN, "r");
   // debug_printf(DP_WARN, "!!!!!!!!!!!%s\n",FN);
	if (fp == NULL) {
        debug_printf(DP_ERROR, "Couldn't open script file!!!");
        exit(EXIT_FAILURE); 
    }

    bool NormalOp=false;

	while ((read = getline(&line, &len, fp)) != -1) {
		if(read>0) {
			if(line[0]=='#') {
				debug_printf(DP_INFO,"%s", line);
                continue;
			}
            
            for(int i = 0; line[i]; i++) { line[i] = tolower(line[i]);  }
			token = strtok(line, " ,.-\n");

            if(token[0]=='_') {  // add
                AddToDict(&token[1], LinopCounter);
                // token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
            }
//             printf( "Token: XX%sXX\n", token );
            if(strcmp(token,"define")==0) {
                char* defname = strtok(NULL, " ,.\n");
				GET_NEXT_TOKEN
                long val=atoi(token);
                debug_printf(DP_INFO,"Defining %s to be %ld\n",defname,val);
                AddToDict(defname, val);
                continue;
			}

            if(strcmp(token,"name")==0) {
                char* defname = strtok(NULL, " ,.\n");
                debug_printf(DP_INFO,"Naming %s to be %ld\n",defname,LinopOutCounter);
                AddToDict(defname, LinopOutCounter);
                continue;
            }
            if(strcmp(token,"setdefop")==0) {
                debug_printf(DP_INFO,"Setting %ld to be the default op\n",LinopOutCounter);
                DefOp=LinopOutCounter;
                continue;
            }

            if(strcmp(token,"fft")==0) {
				//token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                long FFTFlags=atoi(token);
                debug_printf(DP_INFO,"Linop %ld: Adding: FFT with flag %ld\n",LinopCounter,FFTFlags);
                
                NEW_OP_POINTER = linop_fft_create(DIMS, CurDims, FFTFlags);
                ADD_OP
			}
            if(strcmp(token,"ifft")==0) {
				// token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                long FFTFlags=atoi(token);
                debug_printf(DP_INFO,"Linop %ld: Adding: IFFT with flag %ld\n",LinopCounter,FFTFlags);
                
                NEW_OP_POINTER = linop_ifft_create(DIMS, CurDims, FFTFlags);
                ADD_OP
			}
            if(strcmp(token,"fftc")==0) {
				//token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                long FFTFlags=atoi(token);
                debug_printf(DP_INFO,"Linop %ld: Adding: FFTC with flag %ld\n",LinopCounter,FFTFlags);
                
                NEW_OP_POINTER = linop_fftc_create(DIMS, CurDims, FFTFlags);
                ADD_OP
			}
            if(strcmp(token,"ifftc")==0) {
				//token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                long FFTFlags=atoi(token);
                debug_printf(DP_INFO,"Linop %ld: Adding: IFFTC with flag %ld\n",LinopCounter,FFTFlags);
                
                NEW_OP_POINTER = linop_ifftc_create(DIMS, CurDims, FFTFlags);
                ADD_OP
			}
            if(strcmp(token,"grad")==0) {
				//token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                long GradFlags=atoi(token);
                debug_printf(DP_INFO,"Linop %ld: Adding: grad with flag %ld\n",LinopCounter,GradFlags);
                
                NEW_OP_POINTER = linop_gradm_create(DIMS, CurDims, GradFlags);
                // debug_printf(DP_INFO,"SOP %d %d \n",linop_domain(Sop)->N,linop_codomain(Sop)->N);
                // debug_printf(DP_INFO,"Cur %d %d \n",linop_domain(LinopsVec[LinopCounter])->N,linop_codomain(LinopsVec[LinopCounter])->N);
                
                ADD_OP
			}
            if(strcmp(token,"finitediff")==0) {
				//token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                long FDFlags=atoi(token);
                token = strtok(NULL, " ,.");
                long bSnip=atoi(token);

                debug_printf(DP_INFO,"Linop %ld: Adding: finitediff with flag %ld snip %ld\n",LinopCounter,FDFlags, bSnip);
                
                NEW_OP_POINTER = linop_finitediff_create(DIMS, CurDims, FDFlags, bSnip);
                ADD_OP
			}
            if(strcmp(token,"zfinitediff")==0) {
				//token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                long zFDFlags=atoi(token);
                //token = strtok(NULL, " ,.");
                GET_NEXT_TOKEN
                long bCircular=atoi(token);

                debug_printf(DP_INFO,"Linop %ld: Adding: finitediff with flag %ld circular %ld\n",LinopCounter,zFDFlags, bCircular);
                
                NEW_OP_POINTER = linop_zfinitediff_create(DIMS, CurDims, zFDFlags, bCircular);
                ADD_OP
			}

            if(strcmp(token,"fmac")==0) {
				// token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                // long FileIdx=atoi(token);
                long FileIdx=LookInDictOrStr(token);
                // token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                // long SquashFlags=atoi(token);
                long SquashFlags=LookInDictOrStr(token);
                
          //      debug_printf(DP_WARN,"Linop %ld: Adding: FMAC with file #%ld squash flag %ld\n",LinopCounter,FileIdx,SquashFlags);
                
                debug_print_dims(DP_DEBUG1, DIMS, CurDims);
                debug_print_dims(DP_DEBUG1, DIMS, Fdims[FileIdx]);
              
                long MergedDims[DIMS];
                md_merge_dims(DIMS, MergedDims, CurDims, Fdims[FileIdx]);
                long NewDims[DIMS];
                md_select_dims(DIMS, ~SquashFlags, NewDims, MergedDims);

                long CurFlags=md_nontriv_dims(DIMS,CurDims);
                long NewFlags=md_nontriv_dims(DIMS,NewDims);
                long TFlags=md_nontriv_dims(DIMS,Fdims[FileIdx]);
                
                // update CurDims:
                // md_copy_dims(DIMS, CurDims, NewDims);
                // md_select_dims(DIMS, ~SquashFlags, CurDims, MergedDims);
                
                // debug_printf(DP_INFO,"Flags: %ld %ld %ld\n",CurFlags,NewFlags,TFlags);
                NEW_OP_POINTER = linop_fmac_create(DIMS, MergedDims, 
                    ~NewFlags, ~CurFlags, ~TFlags, dataFiles[FileIdx]);
                ADD_OP
			}
         

            if(strcmp(token,"fmaz")==0) {
                GET_NEXT_TOKEN
                long FileIdx=LookInDictOrStr(token);
                GET_NEXT_TOKEN
                long SquashFlags=LookInDictOrStr(token);
                
                debug_printf(DP_INFO,"Linop %ld: Adding: FMAC on CPU with file #%ld squash flag %ld\n",LinopCounter,FileIdx,SquashFlags);
                
                debug_print_dims(DP_INFO, DIMS, CurDims);
                debug_print_dims(DP_INFO, DIMS, Fdims[FileIdx]);
                
                long MergedDims[DIMS];
                md_merge_dims(DIMS, MergedDims, CurDims, Fdims[FileIdx]);
                long NewDims[DIMS];
                md_select_dims(DIMS, ~SquashFlags, NewDims, MergedDims);
                
                long CurFlags=md_nontriv_dims(DIMS,CurDims);
                long NewFlags=md_nontriv_dims(DIMS,NewDims);
                long TFlags=md_nontriv_dims(DIMS,Fdims[FileIdx]);
                
                NEW_OP_POINTER = linop_fmacOnCPU_create(DIMS, MergedDims, 
                    ~NewFlags, ~CurFlags, ~TFlags, dataFiles[FileIdx]);
                ADD_OP
            }
            if(strcmp(token,"transpose")==0) {
				//token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                long dim1=atoi(token);
                //token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                long dim2=atoi(token);
                debug_printf(DP_INFO,"Linop %ld: Adding: Transpose with dims %ld,%ld\n",LinopCounter,dim1,dim2);
                
               // NEW_OP_POINTER = linop_transpose_create(DIMS, CurDims, dim1,dim2);
                NEW_OP_POINTER = linop_transpose_create(DIMS, dim1, dim2, CurDims);
                ADD_OP
			}
            if(strcmp(token,"print")==0) {
				//token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                long msgId=atol(token);
                debug_printf(DP_INFO,"Linop %ld: Adding: print with messageId %ld\n",LinopCounter,msgId);
                
                NEW_OP_POINTER = linop_print_create(DIMS, CurDims, msgId);
                ADD_OP
			}
            if(strcmp(token,"ident")==0) {
				debug_printf(DP_INFO,"Linop %ld: Adding: identity: do nothing\n",LinopCounter);
                
                NEW_OP_POINTER = linop_identity_create(DIMS, CurDims);
                ADD_OP
			}
            if(strcmp(token,"samp")==0) { 
				//token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                long FileIdx=atoi(token);
                
                debug_printf(DP_INFO,"Linop %ld: Adding: Sampling with file #%ld\n",LinopCounter,FileIdx);
                NEW_OP_POINTER = linop_samplingGeneral_create(CurDims, Fdims[FileIdx], dataFiles[FileIdx]);
                ADD_OP
			}
       /*     if(strcmp(token,"part")==0) { 
				token = strtok(NULL, " ,.-");
                long dim1=atoi(token);
                token = strtok(NULL, " ,.-");
                long dim2=atoi(token);
                token = strtok(NULL, " ,.-");
                long K=atoi(token);
                debug_printf(DP_INFO,"Linop %ld: Adding: Partition with dims %ld,%ld and K=%ld\n",LinopCounter,dim1,dim2,K);
                
                NEW_OP_POINTER = linop_PartitionDim_create(DIMS, CurDims, dim1, dim2, K);
                ADD_OP
			} */
       /*     if(strcmp(token,"hankel")==0) { 
				token = strtok(NULL, " ,.-");
                long dim1=atoi(token);
                token = strtok(NULL, " ,.-");
                long dim2=atoi(token);
                token = strtok(NULL, " ,.-");
                long K=atoi(token);
                debug_printf(DP_INFO,"Linop %ld: Adding: Hankel with dims %ld,%ld and K=%ld\n",LinopCounter,dim1,dim2,K);
                
                NEW_OP_POINTER = linop_Hankel_create(DIMS, CurDims, dim1, dim2, K);
                ADD_OP
			} */
            if(strcmp(token,"nufft")==0) {
                //token = strtok(NULL, " ,.");
                GET_NEXT_TOKEN
                long TrajFileIdx=atoi(token);
                //token = strtok(NULL, " ,.");
                GET_NEXT_TOKEN
                long WeightsFileIdx=atoi(token);
                //token = strtok(NULL, " ,.");
                GET_NEXT_TOKEN
                long BasisFileIdx=atoi(token);
                //token = strtok(NULL, " ,.");
                GET_NEXT_TOKEN
                long NUFlags=atoi(token);
                //token = strtok(NULL, " ,.");
                GET_NEXT_TOKEN
                long Toep=atoi(token);
                //token = strtok(NULL, " ,.");
                GET_NEXT_TOKEN
                long pcycle=atoi(token);
                //token = strtok(NULL, " ,.");
                GET_NEXT_TOKEN
                long periodic=atoi(token);
                //token = strtok(NULL, " ,.");
                GET_NEXT_TOKEN
                long lowmem=atoi(token);

                debug_printf(DP_INFO,"Linop %ld: Adding: NUFFT  with file #%ld Weights file %ld Basis file %ld Flags %ld Toeplitz %ld pcycle %ld periodic %ld lowmem %ld\n",
                    LinopCounter,TrajFileIdx,WeightsFileIdx,BasisFileIdx,NUFlags,Toep,pcycle,periodic,lowmem);
                debug_printf(DP_INFO,"----------- NUFFT Trajectory should be [3, Readout, Spokes] !!!\n");
                
                debug_print_dims(DP_INFO, DIMS, CurDims);
                debug_print_dims(DP_INFO, DIMS, Fdims[TrajFileIdx]);
                
                struct nufft_conf_s nuconf = nufft_conf_defaults;
                nuconf.flags = NUFlags;
                nuconf.toeplitz = Toep>0;
                nuconf.pcycle = pcycle>0;
                nuconf.periodic = periodic>0;
                nuconf.lowmem = lowmem>0;

                long ksp_dims[DIMS];
                md_select_dims(DIMS, PHS1_FLAG|PHS2_FLAG, ksp_dims, Fdims[TrajFileIdx]);
                md_copy_dims(DIMS - 3, ksp_dims + 3, CurDims + 3);

                // NEW_OP_POINTER;
                if(BasisFileIdx<0 && WeightsFileIdx<0) {
                    NEW_OP_POINTER= nufft_create(DIMS, ksp_dims, CurDims, Fdims[TrajFileIdx], dataFiles[TrajFileIdx], NULL, nuconf);
                } else if(BasisFileIdx<0 && WeightsFileIdx>=0) {
                    NEW_OP_POINTER= nufft_create2(DIMS, ksp_dims, CurDims, Fdims[TrajFileIdx], dataFiles[TrajFileIdx], Fdims[WeightsFileIdx],dataFiles[WeightsFileIdx],NULL,NULL, nuconf);
                } else if(BasisFileIdx>=0 && WeightsFileIdx<0) {
                    NEW_OP_POINTER= nufft_create2(DIMS, ksp_dims, CurDims, Fdims[TrajFileIdx], dataFiles[TrajFileIdx], NULL,NULL, Fdims[BasisFileIdx],dataFiles[BasisFileIdx], nuconf);
                } else if(BasisFileIdx>=0 && WeightsFileIdx>=0) {
                    NEW_OP_POINTER= nufft_create2(DIMS, ksp_dims, CurDims, Fdims[TrajFileIdx], dataFiles[TrajFileIdx],
                                        Fdims[WeightsFileIdx],dataFiles[WeightsFileIdx], Fdims[BasisFileIdx],dataFiles[BasisFileIdx], nuconf);
                }
                ADD_OP
            
                GET_NEXT_TOKEN
                if(token) { 
                    debug_printf(DP_INFO,"TOKEN %s TOKEN\n",token);
                 //   save_nufft_toep_data(LinopsVec[LinopCounter-1],token);
                } else {
                    debug_printf(DP_INFO,"NUFFT: No other token\n");
                }
                token=line;
            }

            if(strcmp(token,"wavelet")==0) { 
				//token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                long xflags=atoi(token);
                //token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                long randshift=atoi(token);
                debug_printf(DP_INFO,"Linop %ld: Adding: Wavelet with xflags %ld randshift %ld\n",LinopCounter,xflags,randshift);

                long cur_strs[DIMS];
			    md_calc_strides(DIMS, cur_strs, CurDims, CFL_SIZE);

                long minsize[DIMS] = { [0 ... DIMS - 1] = 1 };
                minsize[0] = MIN(CurDims[0], 16);
                minsize[1] = MIN(CurDims[1], 16);
                minsize[2] = MIN(CurDims[2], 16);


                unsigned int wflags = 0;
                unsigned int wxdim = 0;
                for (unsigned int i = 0; i < DIMS; i++) {

                    if ((1 < CurDims[i]) && MD_IS_SET(xflags, i)) {

                        wflags = MD_SET(wflags, i);
                        minsize[i] = MIN(CurDims[i], 16);
                        wxdim += 1;
                    }
                }
                enum wtype wtype = WAVELET_DAU2;
                NEW_OP_POINTER = linop_wavelet_create(DIMS, wflags, CurDims, cur_strs, wtype, minsize, randshift>0);
                ADD_OP
			}
            if(strcmp(token,"resize")==0) { 
				//token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                long FileIdx=atoi(token);
                
                debug_printf(DP_INFO,"Linop %ld: Adding: resize by #%ld\n",LinopCounter,FileIdx);
                NEW_OP_POINTER = linop_resize_create(DIMS, dataFiles[FileIdx],CurDims);
                ADD_OP
			}
            if(strcmp(token,"dblszc")==0) { 
				//token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                long dblFlags=atoi(token);

                long dbl_dims[DIMS];
                md_copy_dims(DIMS,dbl_dims,CurDims);
                for (unsigned int i = 0; i < DIMS; i++) {
                    if (MD_IS_SET(dblFlags, i)) {
                        dbl_dims[i] *= 2;
                    }
                }

                debug_printf(DP_INFO,"Linop %ld: Adding: dblszc flags  #%ld\n",LinopCounter,dblFlags);

                debug_printf(DP_INFO,"dbl_dims: ");
                debug_print_dims(DP_INFO,DIMS,dbl_dims);
                
                NEW_OP_POINTER = linop_resize_create(DIMS, dbl_dims,CurDims);
                ADD_OP
			}
            if(strcmp(token,"halfszc")==0) { 
				//token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                long halfFlags=LookInDictOrStr(token);

                long half_dims[DIMS];
                md_copy_dims(DIMS,half_dims,CurDims);
                for (unsigned int i = 0; i < DIMS; i++) {
                    if (MD_IS_SET(halfFlags, i)) {
                        half_dims[i] /= 2;
                    }
                }

                debug_printf(DP_INFO,"Linop %ld: Adding: halfszc flags  #%ld\n",LinopCounter,halfFlags);

                debug_printf(DP_INFO,"half_dims: ");
                debug_print_dims(DP_INFO,DIMS,half_dims);
                
                NEW_OP_POINTER = linop_resize_create(DIMS, half_dims,CurDims);
                ADD_OP
			}
            if(strcmp(token,"dblsz")==0) { 
				//token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                long dblFlags=atoi(token);

                long dbl_dims[DIMS];
                md_copy_dims(DIMS,dbl_dims,CurDims);
                for (unsigned int i = 0; i < DIMS; i++) {
                    if (MD_IS_SET(dblFlags, i)) {
                        dbl_dims[i] *= 2;
                    }
                }

                debug_printf(DP_INFO,"Linop %ld: Adding: dblsz flags  #%ld\n",LinopCounter,dblFlags);

                debug_printf(DP_INFO,"dbl_dims: ");
                debug_print_dims(DP_INFO,DIMS,dbl_dims);
                
                NEW_OP_POINTER = linop_resizeBase_create(DIMS, dbl_dims,CurDims);
                ADD_OP
			}
            if(strcmp(token,"halfsz")==0) { 
				//token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                long halfFlags=LookInDictOrStr(token);

                long half_dims[DIMS];
                md_copy_dims(DIMS,half_dims,CurDims);
                for (unsigned int i = 0; i < DIMS; i++) {
                    if (MD_IS_SET(halfFlags, i)) {
                        half_dims[i] /= 2;
                    }
                }

                debug_printf(DP_INFO,"Linop %ld: Adding: halfsz flags  #%ld\n",LinopCounter,halfFlags);

                debug_printf(DP_INFO,"half_dims: ");
                debug_print_dims(DP_INFO,DIMS,half_dims);
                
                NEW_OP_POINTER = linop_resizeBase_create(DIMS, half_dims,CurDims);
                ADD_OP
			}

            if(strcmp(token,"byslice")==0) { 
                GET_NEXT_TOKEN
                long WhichOp=LookInDictOrStr(token);
                GET_NEXT_TOKEN
                long WhichDim=LookInDictOrStr(token);

                GET_NEXT_TOKEN
                // debug_printf(DP_INFO,"Token XX%sXX\n",token);
                
                complex float* dataFilesx[MAX_TENSORS_BYSLICE];
                long tstrsx[MAX_TENSORS_BYSLICE];
                const struct linop_s* Opsx[MAX_TENSORS_BYSLICE];

                long tstrs[DIMS];

                // debug_printf(DP_INFO,"XXX\n");

                long fOp, fFile;
                long nTnsrs=0;
                while( token != NULL ) {
                    debug_printf(DP_INFO,"Token XX%sXX\n",token);
                    // debug_printf(DP_INFO,"AAA %ld\n",nTnsrs);
                    fOp=LookInDictOrStr(token);
                    GET_NEXT_TOKEN
                    fFile=LookInDictOrStr(token);
                    debug_printf(DP_INFO,"\t\tAlso op/file %ld/%ld\n",fOp,fFile);
                    
                    Opsx[nTnsrs]=LinopsVec[fOp];
                    dataFilesx[nTnsrs]=dataFiles[fFile];
                    md_calc_strides(DIMS, tstrs, Fdims[fFile], CFL_SIZE);
                    tstrsx[nTnsrs]=tstrs[WhichDim]/CFL_SIZE;
                    
                    GET_NEXT_TOKEN
                    nTnsrs++;
                }

                debug_printf(DP_INFO,"Linop %ld: Adding: byslice OpScript %ld on dim %ld. #Tensors: %ld\n",LinopCounter,WhichOp,WhichDim,nTnsrs);

                
                NEW_OP_POINTER = linop_applyBySlice_create(DIMS, CurDims, WhichDim, LinopsOutVec[WhichOp],dataFilesx,tstrsx,Opsx,nTnsrs);
                // NEW_OP_POINTER = linop_applyBySlice_create(DIMS, CurDims, WhichDim, LinopsOutVec[WhichOp],NULL,NULL,NULL,nTnsrs);
                ADD_OP
                continue;
            }

            if(strcmp(token,"f")==0) {
                //token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
              //  long Idx=atoi(token);
                long Idx=LookInDictOrStr(token);

                debug_printf(DP_INFO,"Adding forward of linop #%ld\n",Idx);

                debug_printf(DP_INFO,"CurDims: ");
                debug_print_dims(DP_INFO,DIMS,CurDims);

                debug_printf(DP_INFO,"LinopIdx dims: ");
                debug_print_dims(DP_INFO,DIMS,linop_domain(LinopsVec[Idx])->dims);
                debug_print_dims(DP_INFO,DIMS,linop_codomain(LinopsVec[Idx])->dims);
                
                Sop = linop_chain(Sop,LinopsVec[Idx]);
                md_copy_dims(DIMS, CurDims, linop_codomain(Sop)->dims);
            }
            if(strcmp(token,"a")==0) {
                //token = strtok(NULL, " ,.-");

                GET_NEXT_TOKEN
                // long Idx=atoi(token);
                long Idx=LookInDictOrStr(token);
              //  debug_printf(DP_WARN,"Adding adjoint of linop #%ld\n",Idx);

                debug_printf(DP_INFO,"CurDims: ");
                debug_print_dims(DP_INFO,DIMS,CurDims);

                debug_printf(DP_INFO,"LinopIdx dims: ");
                debug_print_dims(DP_INFO,DIMS,linop_domain(LinopsVec[Idx])->dims);
                debug_print_dims(DP_INFO,DIMS,linop_codomain(LinopsVec[Idx])->dims);

                Sop = linop_chainAdj(Sop,LinopsVec[Idx]);
                md_copy_dims(DIMS, CurDims, linop_codomain(Sop)->dims);
             //   debug_printf(DP_WARN, "!!!!!!!!!!! \n");

            }
            if(strcmp(token,"n")==0) {
                //token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                long Idx=atoi(token);
                debug_printf(DP_INFO,"Adding Normal of linop #%ld\n",Idx);
                Sop = linop_chainNormal(Sop,LinopsVec[Idx]);
                md_copy_dims(DIMS, CurDims, linop_codomain(Sop)->dims);
            }

            if(strcmp(token,"fx")==0) {
                // token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                //long Idx=atoi(token);
                long Idx=LookInDictOrStr(token);

                debug_printf(DP_INFO,"Adding forward of linopScript #%ld\n",Idx);

                debug_printf(DP_INFO,"CurDims: ");
                debug_print_dims(DP_INFO,DIMS,CurDims);

                debug_printf(DP_INFO,"LinopIdx dims: ");
                debug_print_dims(DP_INFO,DIMS,linop_domain(LinopsOutVec[Idx])->dims);
                debug_print_dims(DP_INFO,DIMS,linop_codomain(LinopsOutVec[Idx])->dims);
                
                Sop = linop_chain(Sop,LinopsOutVec[Idx]);
                md_copy_dims(DIMS, CurDims, linop_codomain(Sop)->dims);
            }
            if(strcmp(token,"ax")==0) {
                // token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                // long Idx=atoi(token);
                long Idx=LookInDictOrStr(token);

                debug_printf(DP_INFO,"Adding adjoint of linopScript #%ld\n",Idx);

                debug_printf(DP_INFO,"CurDims: ");
                debug_print_dims(DP_INFO,DIMS,CurDims);

                debug_printf(DP_INFO,"LinopIdx dims: ");
                debug_print_dims(DP_INFO,DIMS,linop_domain(LinopsOutVec[Idx])->dims);
                debug_print_dims(DP_INFO,DIMS,linop_codomain(LinopsOutVec[Idx])->dims);

                Sop = linop_chainAdj(Sop,LinopsOutVec[Idx]);
                md_copy_dims(DIMS, CurDims, linop_codomain(Sop)->dims);
            }
            if(strcmp(token,"nx")==0) {
                // token = strtok(NULL, " ,.-");
                GET_NEXT_TOKEN
                long Idx=atoi(token);
                debug_printf(DP_INFO,"Adding Normal of linopScript #%ld\n",Idx);
                Sop = linop_chainNormal(Sop,LinopsOutVec[Idx]);
                md_copy_dims(DIMS, CurDims, linop_codomain(Sop)->dims);
            }

            if(strcmp(token,"normal")==0) {
                debug_printf(DP_INFO,"---------\nMoving to normal\n----------------\n");
                NormalOp=true;
                SopBackup=Sop;
                md_copy_dims(DIMS, CurDims, StartDims[LinopOutCounter]);
                Sop = linop_identity_create(DIMS, CurDims);
            }
            if(strcmp(token,"nextlinop")==0) {
               // debug_printf(DP_WARN,"---------\nMoving to next linop\n----------------\n");
                if(NormalOp) {
                    const struct linop_s* tmp=Sop;
                    Sop=linop_PutFowrardOfBInNormalOfA(SopBackup, tmp);
                    linop_free(tmp);
                    NormalOp=false;
                }
                LinopsOutVec[LinopOutCounter]=Sop;
                LinopOutCounter++;
                md_copy_dims(DIMS, CurDims, StartDims[LinopOutCounter]);
                Sop = linop_identity_create(DIMS, CurDims);
            }
		}
	}
	fclose(fp);
	if (line) {
		free(line); }

    if(NormalOp) {
        const struct linop_s* tmp=Sop;
        Sop=linop_PutFowrardOfBInNormalOfA(SopBackup, tmp);
        linop_free(tmp);
    }

    LinopsOutVec[LinopOutCounter]=Sop;


    // FreeLinops();
    debug_printf(DP_DEBUG1,"getLinopScriptFromFile end\n");
    return LinopsOutVec[DefOp];
}

static const char usage_str[] = "<OpScriptTxt> <StartDims> <input> [<file0> [<file> [<file2> [...]]]] <output>";
static const char help_str[] =
		"Apply linop from script -\n"
        "linopScript <OpScriptTxt> <StartDims> <input> [<file0> [<file> [<file2> [...]]]] <output>\n"
		"-----------------------------------------\n"
		"Apply operator script from OpScriptTxt on the input, and save in output\n"
        "Uses other files if mentioned\n"
        "Linops:\n"
        "FFT/IFFT/FFTC/IFFTC <FFT_FLAGS>\n"
        "FMAC <Which_file_no> <SQUASH_FLAGS> : multiplies and then sums\n"
        "Transpose <dim1> <dim2> : transposes the dims\n"
        "Print <messageId> : print messageId on frwrd/adjoint/normal calls\n"
        "ident - do nothing\n"
        "Samp <Which_file_no> : Sampling is multiplication by binary map - so forward=adjoint=normal\n"
        "Part Dim1 Dim2 K\n"
        "Hankel Dim1 Dim2 K\n"
        "resize FileIdx : File contains new size -- not tested\n"
        "dblsz/halfsz/dblszc/halfszc Flags : Flags for dims to double/half\n"
        "wavelet xflags randshift : randshift BOOL 0/1\n"
        "NUFFT TrajFileIdx WeightsFileIdx BasisFileIdx NUFlags ToepBool pcycleBool periodicBool lowmemBool : TrajFile should be [3 readout spokes]. Bool 0/1. NUFT defaults are false.\n"
        "GRAD <Grad_flags> : *Note*: Grad results goes to the 16th dim, not BART's usual 17th\n"
        "finitediff <FD_flags> <boolean snip 0/1>\n"
        "zfinitediff <zFD_flags> <boolean circular 0/1>\n"
        "f/a/n OpNumber\n"
        "fx/ax/nx OpScriptNumber\n"
        "byslice OpScriptNumber WhichDim\n"
        "define string number :  let this string mean this number\n"
        "name string : give name to the current OpScript\n"
        "SetDefOp : set this OpScript as the default one";

int main_linopScript(int argc, char* argv[])
{
    LinopCounter=0;

	bool Normal = false;
    bool Adj = false;
    unsigned int fftmod_flags = 0;

    unsigned int WhichLinopToRun = 0;
    unsigned int nIters = 1;

    bool gpu = false;
    unsigned int gpun = 0;

	const struct opt_s opts[] = {
		OPT_SET('N', &Normal, "Apply normal"),
        OPT_SET('A', &Adj, "Apply adjoint"),
        //OPT_UINT('j', &fftmod_flags, "fftmod_flags", "flags for fftmod_flags of (k-space?) input (if forward/normal) and (sensitivity?) file0"),
        OPT_UINT('L', &WhichLinopToRun, "idx", "Which Linop to run"),
        OPT_UINT('n', &nIters, "iters", "# times to repeat the linop (For timing tests)"),
        OPT_INT('d', &debug_level, "level", "Debug level"),
        OPT_SET('g', &gpu, "use GPU"),
        OPT_UINT('G', &gpun, "gpun", "use GPU device gpun"),
	};	

	num_init();
    
   // cmdline(&argc, argv, 0, 100, usage_str, help_str, ARRAY_SIZE(opts), opts);
    cmdline(&argc, argv, ARRAY_SIZE(usage_str), usage_str, help_str, ARRAY_SIZE(opts), opts);

	int num_args = argc - 1;
    
    if (gpu)
        num_init_gpu_device(gpun);
    else
        num_init();
    debug_printf(DP_INFO,"main_linopScript\n");
    debug_printf(DP_INFO,"AAA %d\n",num_args);
    debug_printf(DP_INFO,"%s\n",argv[0]); // "linopScript"
    debug_printf(DP_INFO,"%s\n",argv[1]); // script file
    debug_printf(DP_INFO,"%s\n",argv[2]); // input dims
    debug_printf(DP_INFO,"%s\n",argv[3]); // input
    debug_printf(DP_INFO,"Out: %s\n",argv[num_args]); // input
    debug_printf(DP_INFO,"BBB\n");

    // #ifdef USE_CUDA
//     print_cuda_meminfo();
// #endif
    if(gpu) {
        ReadScriptFiles_gpu(&argv[4],num_args-4);
    } else {
        ReadScriptFiles(&argv[4],num_args-4);
    }

    
    // fftmod(DIMS, getFdims(1), fftmod_flags, getDataFile(1), getDataFile(1));
    
    debug_printf(DP_INFO,"input dims: %s\n",argv[3]);
    debug_printf(DP_INFO,"XXXXXXXXXXXXXX\n");
    long inputDims_dims[DIMS];
    complex float* inputDimsC = load_cfl(argv[2], DIMS, inputDims_dims);

    
    debug_printf(DP_INFO,"inputDims_dims: ");
    debug_print_dims(DP_INFO,DIMS,inputDims_dims);
    
    if(inputDims_dims[0]!=DIMS) { 
        debug_printf(DP_ERROR,"inputDims_dims[0]!=%d\n",DIMS);
        exit(0);
    }
    
    debug_printf(DP_DEBUG1,"Reading script:\n");
    getLinopScriptFromFile(argv[1],inputDimsC,inputDims_dims[1]);
    const struct linop_s* Sop = getLinopOutX(WhichLinopToRun);
    
    long CurDims[DIMS];
    debug_printf(DP_DEBUG1,"WhichLinopToRun %d\n",WhichLinopToRun);
	md_copy_dims(DIMS, CurDims, linop_domain(Sop)->dims);
    debug_printf(DP_DEBUG1,"CurDims: ");
    debug_print_dims(DP_DEBUG1,DIMS,CurDims);
            
    long inputF_dims[DIMS];
    complex float* input=load_cfl(argv[3], DIMS, inputF_dims);

    debug_printf(DP_INFO,"inputF_dims: ");
    debug_print_dims(DP_INFO,DIMS,inputF_dims);

    
    /* Ops:
     * v FFT FFT_FLAGS
     * v IFFT FFT_FLAGS
     * v FMAC SQUASH_FLAGS             : md_merge_dims(N, dims, dims1, dims2); md_select_dims(N, ~squash, dimso, dims);
     * v identity
     * v Transpose                       ? void md_copy2(unsigned int D, const long dim[D], const long ostr[D], void* optr, const long istr[D], const void* iptr, size_t size)
     * v sampling
     * v nufft
     * v resize : pad and crop
     * v GRAD
     * v finitediff
     * v zfinitediff
     * 
     * Sum                          : fmax with singleton. Needs flags to squash
     * Permute
     * matrix
     * realval
     * wavelet
     * cdiag static struct linop_s* linop_gdiag_create(unsigned int N, const long dims[N], unsigned int flags, const complex float* diag, bool rdiag)
     * rdiag
     * conv
     * cdf97
     * nudft
     */
        
    complex float* outa;
    complex float* out;
    
    const struct operator_s* OpToApply;
    long *InDims;
    long *OutDims;

    long CoDims[DIMS];
    md_copy_dims(DIMS, CoDims, linop_codomain(Sop)->dims);        

    if(!Adj && !Normal)  // Forward: Imag->k
    {
        OpToApply=Sop->forward;
        InDims=CurDims;
        OutDims=CoDims;
    } else if(!Normal) { // adjoint
        OpToApply=Sop->adjoint;
        InDims=CoDims;
        OutDims=CurDims;
    } else { // Normal
        OpToApply=Sop->normal;
        InDims=CurDims;
        OutDims=CurDims;
    }

    outa = create_cfl(argv[num_args], DIMS, OutDims);
    if(gpu) {
#ifdef USE_CUDA
        out=md_gpu_move(DIMS, OutDims, outa, CFL_SIZE);
#endif
    } else {
        out=outa; }

// #ifdef USE_CUDA
//     print_cuda_meminfo();
// #endif

    debug_printf(DP_DEBUG1,"Applying the operator !\n");
    for(long i=0;i<nIters;i++) {
// #ifdef USE_CUDA
//         print_cuda_meminfo();
// #endif
        operator_apply(OpToApply,DIMS, OutDims, out, DIMS, InDims, input); }

// #ifdef USE_CUDA
//     print_cuda_meminfo();
// #endif

    if(gpu) {
        md_copy(DIMS, OutDims, outa, out, CFL_SIZE);
    }
   
    debug_printf(DP_INFO,"Saving output\n");
    if(!Adj) {
        unmap_cfl(DIMS, linop_codomain(Sop)->dims, out);
    } else {
        unmap_cfl(DIMS, linop_domain(Sop)->dims, out);
    }
    
    debug_printf(DP_INFO,"inputF_dims: ");
    debug_print_dims(DP_INFO,DIMS,inputF_dims);

    unmap_cfl(DIMS, inputF_dims, input);
    
    ClearReadScriptFiles(&argv[4],num_args-4);
    FreeLinops(false);
    FreeOutLinops();
    
	exit(EXIT_SUCCESS);

	return 0;
}

