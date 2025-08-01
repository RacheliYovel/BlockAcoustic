using LinearAlgebra
using SparseArrays
using jInv.Mesh
using Helmholtz
using Helmholtz.ElasticHelmholtz
using KrylovMethods
using Multigrid

include("getModels.jl")

const TYPE = ComplexF64;
const ITYPE = Int64;


function getAndSolveElasticHelmholtzMonolithic(n,model,shift,omega; levels = 2, inner = 5)

    rho,lambda,mu,M = getModel(model, n);

    ############### get operator and RHS ##################        

    aten = 0.01*pi; # physical attenuation
    alpha = aten + shift*omega;
    pad = 20; # padding for absorbing BC

    neumanOnTop = true;
    gamma = getABL(M.n, neumanOnTop, [pad; pad], omega) .+ aten;
    gamma_shifted = getABL(M.n, neumanOnTop, [pad; pad], omega) .+ alpha;
    param = ElasticHelmholtzParam(M, omega, lambda, rho, mu, gamma, neumanOnTop, true)
    param_shifted = ElasticHelmholtzParam(M, omega, lambda, rho, mu, gamma_shifted, neumanOnTop, true)
    ElasticOperator = GetElasticHelmholtzOperator(param);
    ElasticOperatorShifted = GetElasticHelmholtzOperator(param_shifted);

    ElasticRHS = getElasticPointSource(M,ComplexF64);
    ElasticRHS = vec(ElasticRHS);
    b = [ElasticRHS; zeros(prod(M.n))];

    ############### setup for the MG ##################

    relaxType = "VankaFaces";
    relaxParam = [0.65;0.5;0.3];
    nu = [1;1];
    println("pre+post relax num is: ",nu[1]+nu[2])

    numCores 	= 4; 
    maxIter     = 200;
    relativeTol = 1e-6;
    cycleType   ='W';
    coarseSolveType = "NoMUMPS";
                
    MG = getMGparam(TYPE, ITYPE, levels,numCores,maxIter,relativeTol,relaxType,relaxParam,nu[1],nu[2],cycleType,coarseSolveType,0.0,0.0,"SystemsFacesMixedLinear");
    MG = MGsetup(ElasticOperatorShifted,M,MG,1,false);
    
    ####################################################
    
    x,param,iter,resvec = solveGMRES_MG(ElasticOperator,MG,b,zeros(eltype(b),size(b)),true,inner,true);

end
