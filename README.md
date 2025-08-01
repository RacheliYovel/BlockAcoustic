# BlockAcoustic
A block-acoustic preconditioner for the elastic Helmholtz equation. 
Rachel Yovel, Eran Treister
contact: yovelr@bgu.ac.il

This repository includes the code for a block triangular preconditioner for elastic Helmholtz whose blocks are acoustic Helmholtz operators, following the paper https://arxiv.org/abs/2411.15897

To reproduce the results from the paper, run the file "driver.jl". To enable the experiments with advanced heterogeneous models, make sure that you have a local folder named "BenchmarkModels" including the SEG/EAGE salt model, SEAM Phase I model, Marmousi2 and Overthrust 3D model, see references in the paper.

This code relies on the following Julia packages:
LinearAlgebra
SparseArrays

And the following JuliaInv packages  (see https://github.com/JuliaInv)
jInv.Mesh
KrylovMethods
Multigrid
Helmholtz
Helmholtz.ElasticHelmholtz

The code is organized as follows:

The file "getBlockPrec.jl" contains functions that produce the elastic Helmholtz problem in mixed formulation and all the blocks of the block-preconditioner, and some setup for the solution of each block within the block preconditioner (the LU decomposition of each block in a case of direct solution, and multigrid setup in a case of multigrid solution). This file elables also a simple finite difference implementation of older block preconditioners: Fp (PCD - pressure convection diffusion preconditioning) and BFBt (LSC - least square commutator preconditioning) of Elman et. al.

The file "getModels.jl" produces density and Lame coefficients for the constant and heterogeneous velocity models.

The file "multigrid.jl" includes a multigrid framework that enables the solution of different staggering (faces, edges, centers etc.), with different intergrid and smoother options, different cycle types (V,W or more recursive calls) and an exact or inexact coarse solve.

The file "solveGMRES_blockprec.jl" builds the block preconditioner's function and solves the system by preconditioned GMRES.

Finally, for the sake of comparison with monolithic multigrid from the paper https://arxiv.org/abs/1806.11277, the file "monolithic.jl" calls code from the jInv Helmholtz and Multigrid packages above to solve the elastic Helmholtz equation in mixed formulation by monolithic multigrd.
