# ssamp
README for demonstration of compressed sensing (CS) recovery with 1D finite-difference sparsity
by Jaewook Kang (jwkkang@gist.ac.kr) 
=============================================================

Description: 
=========================================================================
 This demonstration contains methods for performing CS recovery with 1D finite-difference
 
 The solvers included in this comparison are as given below:

	-ssAMP-BGFD : Jaewook Kang, Hyoyoung Jung, Heung-No Lee, and Kiseon Kim,  
               "Bernoulli-Gaussian Approximate Message-Passing Algorithm for Compressed 
                Sensing with 1D-Finite-Difference Sparsity,"  submitted SEP 2015
               (https://sites.google.com/site/jwkang10/)
               
	-EFLA :  J. Liu, L. Yuan, and J. Ye. 
                 An efficient algorithm for a class of fused lasso problems,
                 proc of ACM SIGKDD Conference on Knowledge Discovery and Data Mining, 2010
                (http://yelab.net/software/SLEP/ )
                
	-TVAMP:  D. L. Donoho, I. Johnstone, and A. Montanari, 
                  Accurate prediction of phase transitions in compressed sensing via 
                  a connection to minimax denoising,  IEEE Trans. Inform. Theory, 
                  vol. 59, no. 6, pp. 3396-3433, June 2013.
                  
	-Chambolle-Pock : A. Chambolle, T. Pock, 
                  A first-order primal-dual algorithm for convex
                  problems with applications to imaging, J. Math. Imag. Vis., vol. 40, pp.
                  120-145, May 2011
                  
	-GrAMPA :  M. Borgerding and P. Schniter, 
                Generalized approximate message passing for the cosparse analysis
                model,ICASSP 2015 and avabilable at ArXiv:1312.3968v1 [cs.IT], Dec. 2013
                ( http://www2.ece.ohio-state.edu/~schniter/GrAMPA/index.html ) 
		
==============================================================================

This demonstration require two external numerical solvers for the TVAMP denoiser

1) condat1DTV.c (MEX version) by L Condat
available at http://www.gipsa-lab.grenoble-inp.fr/~laurent.condat/publications.html 

2) flsa.h and flsa.c by J. Liu
available at http://yelab.net/software/SLEP/ 

===============================================================================

Player Information:
=====================
"	MATLAB Version: 8.2.0.701 (R2013b)
"	Operating System: Microsoft Windows 7 Version 6.1 (Build 7601: Service Pack 1)
"	Java Version: Java 1.7.0_11-b21 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode

Packing List: 
  *Demonstration
 	-demo_1DPWC.m          : providing a  comparison among recent algorithms for CS recovery with 1D finite-difference
 	-demo_RIPmatrices.m    : providing a comparison of the ssAMP-BGFD recovery with several measurement matrix satisfying the standard RIP requirement   
 	-demo_CPUruntime.m     : providing a CPU runtime comparison among recent algorithms  for CS recovery with 1D finite-difference
 	-demo_NMSE_over_iter.m : providing a  NMSE over iteration comparison among recent algorithms for CS recovery with 1D finite-difference
 	-demo_SNP.m            : an exemplary CS recovery to SNP genomic data
 	
 	* Drawing Phase transition curves	
 	-PTcurve_ssAMP.m
 	-PTcurve_TVAMP.m
 	-PTcurve_TVCP.m
 	-PTcurve_GrAMPABG.m
 	-PTcurve_EFLA.m
 	
 	*Solvers 
 	- solve_ssAMP.m
 	- solve_TVAMP_Condat.m
 	- solve_TVAMP_FLSA.m
 	- solve_chambolle_pock_TV.m
 	- fusedLeastR.m (external solver for EFLA)
 	
 	*EMtuning for ssAMP
 	-init_prior_parameters.m
 	-ssAMP_EMupdate.m
 	
 	*Sythetic measurement generation
 	-sig_1DFD_gen.m
 	-mtxgen.m
 	
 	* solver configulation
 	- EFLA_config.m
 	- config_GrAMPABG.m

 	* The external codes
  - GAMP_codes  ( http://www2.ece.ohio-state.edu/~schniter/GrAMPA/index.html, gamplatlab20141001.zip)
  - condat1DTV (http://www.gipsa-lab.grenoble-inp.fr/~laurent.condat/publications.html )
  - flsa       (http://yelab.net/software/SLEP/)
  

Contact and feedback Information:
=================================
"	Office phone: +82-62-715-2264
"	E-mail      :  jwkkang@gist.ac.kr, jwkang10@gmail.com
%===============================================================================
Copyright (c) 2015, Gwangju Institute of Science and Technology.
All rights reserved.

Contributors: Jaewook Kang, Hyoyoung Jung, Heung-No Lee and Kiseon Kim

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are
met:

	1. Redistributions of source code must retain the above copyright
	notice, this list of conditions and the following disclaimer.

	2. Redistributions in binary form must reproduce the above copyright
	notice, this list of conditions and the following disclaimer in the
	documentation and/or other materials provided with the distribution.

	3. Neither the names of  Gwangju Institute of Science and Technology, 
	nor the names of its contributors may
	be used to endorse or promote products derived from this software
	without specific prior written permission.

===================================================================================




      
