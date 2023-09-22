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
#include "linops/opScript.h"
#include "linops/fmac.h"
#include "linops/sampling.h"
#include "linops/someops.h"
#include "linops/sampling.h"
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

#define MAX_OPS 20
#define MAX_OpFILES 10
#define Op_ADD_LinOP Sop = operator_chain(Sop,NewOp->forward);md_copy_dims(DIMS, CurDims, operator_codomain(Sop)->dims);OpLinopsVec[opLinopCounter++]=NewOp;
// #define ADD_OP const struct linop_s* tmp = Sop;Sop = linop_chain(tmp,NewOp);linop_free(tmp);linop_free(NewOp);md_copy_dims(DIMS, CurDims, linop_codomain(Sop)->dims);

complex float* OpdataFiles[MAX_OpFILES];
long OpFdims[MAX_OpFILES][DIMS];

const struct linop_s* OpLinopsVec[MAX_OPS];
long opLinopCounter=0;

long * getOpFdims(long i) { return OpFdims[i]; }
complex float* getOpDataFile(long i) { return OpdataFiles[i]; }

void OpFreeLinops() {
    debug_printf(DP_INFO,"Freeing %ld Ops\n",opLinopCounter);
    for(long i=0;i<opLinopCounter;i++) {
        linop_free(OpLinopsVec[i]);
    }
    debug_printf(DP_INFO,"FreeOps done\n");
}
void OpReadScriptFiles(char* argv[],long n) {
    long i;
    debug_printf(DP_INFO,"Reading files\n");
    for(i=0;i<n;i++) {        
        OpdataFiles[i] = load_cfl(argv[i], DIMS, OpFdims[i]);
        printf("Reading %s: %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld\n",argv[i],OpFdims[i][0],OpFdims[i][1],OpFdims[i][2],OpFdims[i][3],OpFdims[i][4],OpFdims[i][5],OpFdims[i][6],OpFdims[i][7],OpFdims[i][8],OpFdims[i][9],OpFdims[i][10],OpFdims[i][11],OpFdims[i][12],OpFdims[i][13],OpFdims[i][14],OpFdims[i][15]);
    }
    debug_printf(DP_INFO,"Finished reading files\n");
}

void OpClearReadScriptFiles( char* argv[],long n) {
    long i;
    debug_printf(DP_INFO,"Clearing files' memory\n");
    for(i=0;i<n;i++) {
        debug_printf(DP_INFO,"Clearing %s\n",argv[i]);
        unmap_cfl(DIMS, OpFdims[i], OpdataFiles[i]);
    }
    debug_printf(DP_INFO,"Finished Clearing files' memory\n");
}

const struct operator_s* getOpScriptFromFile(const char *FN, long CurDims[]) {
    opLinopCounter=0;

    const struct operator_s* Sop = operator_identity_create(DIMS, CurDims);

    FILE * fp;
	char * line = NULL;
    char *token;
	size_t len = 0;
	ssize_t read;
    debug_printf(DP_INFO,"getOpScriptFromFile start\n");
    
// 	fp = fopen(argv[1], "r");
    fp = fopen(FN, "r");
	if (fp == NULL) {
        debug_printf(DP_ERROR, "Couldn't open script file!!!");
        exit(EXIT_FAILURE); }
    // debug_printf(DP_INFO,"--------------\n");
    
	while ((read = getline(&line, &len, fp)) != -1) {
// 		debug_printf(DP_INFO,"Retrieved line of length %zu:\n", read);
 		// debug_printf(DP_INFO,"LINE READ: %s", line);
		if(read>0) {
			if(line[0]=='#') {
				debug_printf(DP_INFO,"%s", line);
                continue;
			}
            
            for(int i = 0; line[i]; i++) { line[i] = tolower(line[i]);  }
			token = strtok(line, " ,.-");
//             printf( "Token: XX%sXX\n", token );
            if(strcmp(token,"fft")==0) {
				token = strtok(NULL, " ,.-");
                long FFTFlags=atoi(token);
                debug_printf(DP_INFO,"Linop %ld: Adding: FFT with flag %ld\n",opLinopCounter,FFTFlags);
                
                const struct linop_s* NewOp = linop_fft_create(DIMS, CurDims, FFTFlags);
                Op_ADD_LinOP
			}
            if(strcmp(token,"ifft")==0) {
				token = strtok(NULL, " ,.-");
                long FFTFlags=atoi(token);
                debug_printf(DP_INFO,"Linop %ld: Adding: IFFT with flag %ld\n",opLinopCounter,FFTFlags);
                
                const struct linop_s* NewOp = linop_ifft_create(DIMS, CurDims, FFTFlags);
                Op_ADD_LinOP
			}
            if(strcmp(token,"fftc")==0) {
				token = strtok(NULL, " ,.-");
                long FFTFlags=atoi(token);
                debug_printf(DP_INFO,"Linop %ld: Adding: FFTC with flag %ld\n",opLinopCounter,FFTFlags);
                
                const struct linop_s* NewOp = linop_fftc_create(DIMS, CurDims, FFTFlags);
                Op_ADD_LinOP
			}
            if(strcmp(token,"ifftc")==0) {
				token = strtok(NULL, " ,.-");
                long FFTFlags=atoi(token);
                debug_printf(DP_INFO,"Linop %ld: Adding: IFFTC with flag %ld\n",opLinopCounter,FFTFlags);
                
                const struct linop_s* NewOp = linop_ifftc_create(DIMS, CurDims, FFTFlags);
                Op_ADD_LinOP
			}
            if(strcmp(token,"fmac")==0) {
				token = strtok(NULL, " ,.-");
                long FileIdx=atoi(token);
                token = strtok(NULL, " ,.-");
                long SquashFlags=atoi(token);
                
                debug_printf(DP_INFO,"Linop %ld: Adding: FMAC with file #%ld squash flag %ld\n",opLinopCounter,FileIdx,SquashFlags);
                
                debug_print_dims(DP_INFO, DIMS, CurDims);
                debug_print_dims(DP_INFO, DIMS, OpFdims[FileIdx]);
                
                long MergedDims[DIMS];
                md_merge_dims(DIMS, MergedDims, CurDims, OpFdims[FileIdx]);
                long NewDims[DIMS];
                md_select_dims(DIMS, ~SquashFlags, NewDims, MergedDims);
                
                long CurFlags=md_nontriv_dims(DIMS,CurDims);
                long NewFlags=md_nontriv_dims(DIMS,NewDims);
                long TFlags=md_nontriv_dims(DIMS,OpFdims[FileIdx]);
                
                // update CurDims:
                md_copy_dims(DIMS, CurDims, NewDims);
                // md_select_dims(DIMS, ~SquashFlags, CurDims, MergedDims);
                
                // debug_printf(DP_INFO,"Flags: %ld %ld %ld\n",CurFlags,NewFlags,TFlags);
                const struct linop_s* NewOp = linop_fmac_create(DIMS, MergedDims, 
                    ~NewFlags, ~CurFlags, ~TFlags, OpdataFiles[FileIdx]);
                Op_ADD_LinOP
			}
            if(strcmp(token,"transpose")==0) {
				token = strtok(NULL, " ,.-");
                long dim1=atoi(token);
                token = strtok(NULL, " ,.-");
                long dim2=atoi(token);
                debug_printf(DP_INFO,"Linop %ld: Adding: Transpose with dims %ld,%ld\n",opLinopCounter,dim1,dim2);
                
                const struct linop_s* NewOp = linop_transpose_create(DIMS, CurDims, dim1,dim2);
                Op_ADD_LinOP
			}
            if(strcmp(token,"print")==0) {
				token = strtok(NULL, " ,.-");
                long msgId=atol(token);
                debug_printf(DP_INFO,"Linop %ld: Adding: print with messageId %ld\n",opLinopCounter,msgId);
                
                const struct linop_s* NewOp = linop_print_create(DIMS, CurDims, msgId);
                Op_ADD_LinOP
			}
            if(strcmp(token,"ident")==0) {
				debug_printf(DP_INFO,"Linop %ld: Adding: identity: do nothing\n",opLinopCounter);
                
                const struct linop_s* NewOp = linop_identity_create(DIMS, CurDims);
                Op_ADD_LinOP
			}
            if(strcmp(token,"samp")==0) { 
				token = strtok(NULL, " ,.-");
                long FileIdx=atoi(token);
                
                debug_printf(DP_INFO,"Linop %ld: Adding: Sampling with file #%ld\n",opLinopCounter,FileIdx);
                const struct linop_s* NewOp = linop_samplingGeneral_create(CurDims, OpFdims[FileIdx], OpdataFiles[FileIdx]);
                Op_ADD_LinOP
			}
            if(strcmp(token,"part")==0) { 
				token = strtok(NULL, " ,.-");
                long dim1=atoi(token);
                token = strtok(NULL, " ,.-");
                long dim2=atoi(token);
                token = strtok(NULL, " ,.-");
                long K=atoi(token);
                debug_printf(DP_INFO,"Linop %ld: Adding: Partition with dims %ld,%ld and K=%ld\n",opLinopCounter,dim1,dim2,K);
                
                const struct linop_s* NewOp = linop_PartitionDim_create(DIMS, CurDims, dim1, dim2, K);
                Op_ADD_LinOP
			}
            if(strcmp(token,"hankel")==0) { 
				token = strtok(NULL, " ,.-");
                long dim1=atoi(token);
                token = strtok(NULL, " ,.-");
                long dim2=atoi(token);
                token = strtok(NULL, " ,.-");
                long K=atoi(token);
                debug_printf(DP_INFO,"Linop %ld: Adding: Hankel with dims %ld,%ld and K=%ld\n",opLinopCounter,dim1,dim2,K);
                
              //  const struct linop_s* NewOp = linop_Hankel_create(DIMS, CurDims, dim1, dim2, K);
              //  Op_ADD_LinOP

                // md_copy_dims(DIMS, CurDims, linop_codomain(Sop)->dims);
			}
            if(strcmp(token,"nufft")==0) {
                token = strtok(NULL, " ,.");
                long TrajFileIdx=atoi(token);
                token = strtok(NULL, " ,.");
                long WeightsFileIdx=atoi(token);
                token = strtok(NULL, " ,.");
                long BasisFileIdx=atoi(token);
                token = strtok(NULL, " ,.");
                long NUFlags=atoi(token);
                token = strtok(NULL, " ,.");
                long Toep=atoi(token);
                token = strtok(NULL, " ,.");
                long pcycle=atoi(token);
                token = strtok(NULL, " ,.");
                long periodic=atoi(token);
                token = strtok(NULL, " ,.");
                long lowmem=atoi(token);

                debug_printf(DP_INFO,"Linop %ld: Adding: NUFFT  with file #%ld Weights file %ld Basis file %ld Flags %ld Toeplitz %ld pcycle %ld periodic %ld lowmem %ld\n",
                    opLinopCounter,TrajFileIdx,WeightsFileIdx,BasisFileIdx,NUFlags,Toep,pcycle,periodic,lowmem);
                debug_printf(DP_INFO,"----------- NUFFT Trajectory should be [3, Readout, Spokes] !!!\n");
                
                debug_print_dims(DP_INFO, DIMS, CurDims);
                debug_print_dims(DP_INFO, DIMS, OpFdims[TrajFileIdx]);
                
                struct nufft_conf_s nuconf = nufft_conf_defaults;
                nuconf.flags = NUFlags;
                nuconf.toeplitz = Toep>0;
                nuconf.pcycle = pcycle>0;
                nuconf.periodic = periodic>0;
                nuconf.lowmem = lowmem>0;

                long ksp_dims[DIMS];
                md_select_dims(DIMS, PHS1_FLAG|PHS2_FLAG, ksp_dims, OpFdims[TrajFileIdx]);
                md_copy_dims(DIMS - 3, ksp_dims + 3, CurDims + 3);

                const struct linop_s* NewOp;
                if(BasisFileIdx<0 && WeightsFileIdx<0) {
                    NewOp= nufft_create(DIMS, ksp_dims, CurDims, OpFdims[TrajFileIdx], OpdataFiles[TrajFileIdx], NULL, nuconf);
                } else if(BasisFileIdx<0 && WeightsFileIdx>=0) {
                    NewOp= nufft_create2(DIMS, ksp_dims, CurDims, OpFdims[TrajFileIdx], OpdataFiles[TrajFileIdx], OpFdims[WeightsFileIdx],OpdataFiles[WeightsFileIdx],NULL,NULL, nuconf);
                } else if(BasisFileIdx>=0 && WeightsFileIdx<0) {
                    NewOp= nufft_create2(DIMS, ksp_dims, CurDims, OpFdims[TrajFileIdx], OpdataFiles[TrajFileIdx], NULL,NULL, OpFdims[BasisFileIdx],OpdataFiles[BasisFileIdx], nuconf);
                } else if(BasisFileIdx>=0 && WeightsFileIdx>=0) {
                    NewOp= nufft_create2(DIMS, ksp_dims, CurDims, OpFdims[TrajFileIdx], OpdataFiles[TrajFileIdx],
                                        OpFdims[WeightsFileIdx],OpdataFiles[WeightsFileIdx], OpFdims[BasisFileIdx],OpdataFiles[BasisFileIdx], nuconf);
                }
                Op_ADD_LinOP
            }
            if(strcmp(token,"f")==0) {
                token = strtok(NULL, " ,.-");
                long Idx=atoi(token);
                debug_printf(DP_INFO,"Adding forward of linop #%ld\n",Idx);
                Sop = operator_chain(Sop,OpLinopsVec[Idx]->forward);
                md_copy_dims(DIMS, CurDims, operator_codomain(Sop)->dims);
            }
            if(strcmp(token,"a")==0) {
                token = strtok(NULL, " ,.-");
                long Idx=atoi(token);
                debug_printf(DP_INFO,"Adding adjoint of linop #%ld\n",Idx);
                Sop = operator_chain(Sop,OpLinopsVec[Idx]->adjoint);
                md_copy_dims(DIMS, CurDims, operator_codomain(Sop)->dims);
            }
            if(strcmp(token,"n")==0) {
                token = strtok(NULL, " ,.-");
                long Idx=atoi(token);
                debug_printf(DP_INFO,"Adding Normal of linop #%ld\n",Idx);
                Sop = operator_chain(Sop,OpLinopsVec[Idx]->normal);
                md_copy_dims(DIMS, CurDims, operator_codomain(Sop)->dims);
            }
		}
	}
	fclose(fp);
	if (line) {
		free(line); }
    OpFreeLinops();
    debug_printf(DP_INFO,"getOpScriptFromFile end\n");
    return Sop;
}

static const char usage_str[] = "<OpScriptTxt> <StartDims> <input> [<file0> [<file> [<file2> [...]]]] <output>";
static const char help_str[] =
		"Apply op from script -\n"
        "opScript <OpScriptTxt> <StartDims> <input> [<file0> [<file> [<file2> [...]]]] <output>\n"
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
        "NUFFT TrajFileIdx WeightsFileIdx BasisFileIdx NUFlags ToepBool pcycleBool periodicBool lowmemBool : TrajFile should be [3 readout spokes]. Bool 0/1. NUFT defaults are false.";

int main_opScript(int argc, char* argv[])
{
    opLinopCounter=0;

	unsigned int fftmod_flags = 0;

	const struct opt_s opts[] = {
	    OPT_UINT('j', &fftmod_flags, "fftmod_flags", "flags for fftmod_flags of (k-space?) input (if forward/normal) and (sensitivity?) file0"),
	};	

	num_init();
    
    cmdline2(&argc, argv, 0, 100, usage_str, help_str, ARRAY_SIZE(opts), opts);
	int num_args = argc - 1;
    
    debug_printf(DP_INFO,"main_opScript\n");
    debug_printf(DP_INFO,"AAA %d\n",num_args);
    debug_printf(DP_INFO,"%s\n",argv[0]); // "linopScript"
    debug_printf(DP_INFO,"%s\n",argv[1]); // script file
    debug_printf(DP_INFO,"%s\n",argv[2]); // input dims
    debug_printf(DP_INFO,"%s\n",argv[3]); // input
    debug_printf(DP_INFO,"Out: %s\n",argv[num_args]); // input
    debug_printf(DP_INFO,"BBB\n");
    
    OpReadScriptFiles(&argv[4],num_args-4);
    
    // fftmod(DIMS, getOpFdims(1), fftmod_flags, getDataFile(1), getDataFile(1));
    
    debug_printf(DP_INFO,"input dims: %s\n",argv[3]);
    debug_printf(DP_INFO,"XXXXXXXXXXXXXX\n",argv[3]);
    long inputDims_dims[DIMS];
    long input_dims[DIMS];
    complex float* inputDims = load_cfl(argv[2], DIMS, inputDims_dims);
    
    debug_printf(DP_INFO,"inputDims_dims: ");
    debug_print_dims(DP_INFO,DIMS,inputDims_dims);
    
    md_copy_dims(DIMS, input_dims, inputDims_dims);
    for(long d=0;d<inputDims_dims[0]*inputDims_dims[1];d++) {
        // debug_printf(DP_INFO,"d %d\n",d);
        input_dims[d]=inputDims[d];
    }
    // md_copy_dims(DIMS, input_dims, inputDims);
    debug_printf(DP_INFO,"input_dims: ");
    debug_print_dims(DP_INFO,DIMS,input_dims);
        
    long inputF_dims[DIMS];
    complex float* input=load_cfl(argv[3], DIMS, inputF_dims);;

    long CurDims[DIMS];
	md_copy_dims(DIMS, CurDims, input_dims);
    
    debug_printf(DP_INFO,"Reading script:\n");
    const struct operator_s* Sop =getOpScriptFromFile(argv[1],CurDims);
    
    /* Ops:
     * FFT FFT_FLAGS
     * IFFT FFT_FLAGS
     * FMAC SQUASH_FLAGS             : md_merge_dims(N, dims, dims1, dims2); md_select_dims(N, ~squash, dimso, dims);
     * identity
     * Transpose                       ? void md_copy2(unsigned int D, const long dim[D], const long ostr[D], void* optr, const long istr[D], const void* iptr, size_t size)
     * sampling
     
     * Sum                          : fmax with singleton. Needs flags to squash
     * Permute
     * matrix
     * GRAD?
     * realval
     * wavelet
     * finitediff
     * zfinitediff
     * cdiag static struct linop_s* linop_gdiag_create(unsigned int N, const long dims[N], unsigned int flags, const complex float* diag, bool rdiag)
     * rdiag
     * resize : pad and crop
     * conv
     * cdf97
     * nufft
        struct nufft_conf_s nuconf = nufft_conf_defaults;
        nuconf.toeplitz = true;
        nuconf.lowmem = true;
		OPT_SET('K', &nuconf.pcycle, "randshift for NUFFT"),
        const struct linop_s* fft_op = nufft_create2(DIMS, ksp_dims2, coilim_dims, traj_dims, traj, wgs_dims, weights, basis_dims, basis, conf);
     * nudft
     */
        
    
	complex float* out;
    
    debug_printf(DP_INFO,"Applying the operator\n");
    debug_printf(DP_INFO,"Expected out dims:");
    debug_print_dims(DP_INFO,DIMS,CurDims);
    debug_printf(DP_INFO, "From linop out:");
    debug_print_dims(DP_INFO,operator_codomain(Sop)->N,operator_codomain(Sop)->dims);
    out = create_cfl(argv[num_args], DIMS, CurDims);
    operator_apply(Sop, DIMS, CurDims, out,	DIMS, input_dims, input);

    // To apply a proximal operator
    /*
    long blkdims[DIMS];
    for(long d=0;d<16;d++) {
        blkdims[d]=1;
    }

    int remove_mean = 0;

    bool randshift=false;
    bool overlapping_blocks=false;
    float lambda=0.01;
    unsigned int xflags=51;
    const struct operator_p_s* prox_op=lrthresh_create(CurDims, randshift, xflags, blkdims, lambda, false, remove_mean, overlapping_blocks);

                    // const struct operator_s* opop=operator_p_upcast(prox_op);

    const struct operator_s* opop=opFromOpps(prox_op,0.3);
    operator_apply(opop, DIMS, CurDims, out, DIMS, CurDims, out);

                    // struct iter_op_p_s a_prox_ops = OPERATOR_P2ITOP(prox_op);

                    // a_prox_ops.fun(a_prox_ops.data, 0.3, out, out);
    */



    // fftmod(DIMS, Adjdims, fftmod_flags, adjFile, adjFile);

    // linop_adjoint(Sop, DIMS, input_dims, outadjFile);
   
//    OpFreeLinops();

    debug_printf(DP_INFO,"Saving output\n");
    unmap_cfl(DIMS, operator_codomain(Sop)->dims, out);
    
    OpClearReadScriptFiles(&argv[4],num_args-4);
    unmap_cfl(DIMS, input_dims, input);
    
    // xfree(adj_file);
    
	exit(EXIT_SUCCESS);

	return 0;
}
