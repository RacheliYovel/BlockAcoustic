using LinearAlgebra
using SparseArrays
using jInv.Mesh
using Helmholtz
using Helmholtz.ElasticHelmholtz


function getElasticHelmholtzBlocks(n, model, shift, omega; lambda_factor = 1.0, beta = 1.0)

    ##################### setting the parameters #####################
    
    rho,lambda,mu,M = getModel(model, n);
    n = M.n;

    aten = 0.01 * pi; # physical attenuation
    alpha = aten + shift * omega;
    pad = 20; # padding for absorbing BC
    neumanOnTop = true # compatible with a point source in the middle of upper face

    ##################### size notation #####################
    if length(n) == 2

        nf1 = prod(n + [1; 0])
        nf2 = prod(n + [0; 1])
        nc = prod(n)
        nf = nf1 + nf2; 
        ns = [[n];[nf];[nf1];[nf2];[nc]];
        pad_vec = [pad; pad]
        spread = true;

    elseif length(n) == 3

        nf1 = prod(n + [1; 0; 0])
        nf2 = prod(n + [0; 1; 0])
        nf3 = prod(n + [0; 0; 1])
        nc = prod(n)
        nf = nf1 + nf2 + nf3; 
        ns = [[n];[nf];[nf1];[nf2];[nf3];[nc]];
        pad_vec = [pad; pad; pad]
        spread = false;

    end

    ##################### operator and blocks #####################

    gamma = getABL(n, neumanOnTop, pad_vec, omega)

    param = ElasticHelmholtzParam(M, omega, lambda_factor .* lambda, rho, mu, gamma .+ aten, neumanOnTop, true)
    ElasticOperator = GetElasticHelmholtzOperator(param; spread, beta); # beta = 1.0 for 2nd order MAC, beta = 2/3 for LFA-tuned discretization
    A = ElasticOperator[1:nf,1:nf];
    Bt = ElasticOperator[1:nf,(nf+1):(nf+nc)];
    B = ElasticOperator[(nf+1):(nf+nc),1:nf];
    C = -ElasticOperator[(nf+1):(nf+nc),(nf+1):(nf+nc)];
    ElasticOperator = 0.0;

    ##################### shifted operator and blocks #####################

    paramShifted = ElasticHelmholtzParam(M, omega, lambda_factor .* lambda, rho, mu, gamma .+ alpha, neumanOnTop, true)
    ElasticOperatorShifted = GetElasticHelmholtzOperator(paramShifted; spread, beta);
    lambda = 0.0;

    A1_s = ElasticOperatorShifted[1:nf1,1:nf1];
    A2_s = ElasticOperatorShifted[nf1+1:nf1+nf2,nf1+1:nf1+nf2];
    AsBlocks = [[A1_s];[A2_s]];
    if length(n) == 3
        A3_s = ElasticOperatorShifted[nf1+nf2+1:nf,nf1+nf2+1:nf];
        AsBlocks = [AsBlocks;[A3_s]];
    end

    Shift = blockdiag(AsBlocks...) - A;
    if beta == 1.0
        Shift = diag(Shift);
        Shift = spdiagm(0 => Shift);
    end
    A = 0.0;

    function ElasticOperatorFunc(x)
       xu = x[1:nf];
       xp = x[nf+1:end];
       up = (blockdiag(AsBlocks...) - Shift) * xu + Bt * xp;
       down = B * xu - C * xp;
       return [up;down]
    end

    ##################### RHS point source #####################

    ElasticRHS = getElasticPointSource(M,ComplexF64);
    ElasticRHS = vec(ElasticRHS);
    b = [ElasticRHS; zeros(prod(M.n))];

    absorbing_s = (1.0 .- 1im * (gamma[:] .+ alpha) ./ omega);

    return ElasticOperatorFunc,AsBlocks,Bt,B,C,ns,mu,rho,absorbing_s,b

end


function getPrecBlocks(A,Bt,B,C,mu,rho,omega,absorbing; version = "BlockAc")

    Ap = B * Bt * spdiagm(mu[:]) .- spdiagm(omega^2 .* rho[:] .* absorbing);

    if version == "BlockAc"
        Hp = B*Bt + Ap*C;
    elseif version == "BFBt"
        Hp = B * A * Bt;
    elseif version == "Fp"
        Hp = Ap;
    end

    function ApFunc(x)
        return (Hp - B*Bt) * spdiagm(0 => (1 ./ diag(C))) * x
    end

    Ap = 0.0;

    return Hp, ApFunc

end

function BlockPrecSetup(AsBlocks,Bt,B,C,mu,rho,omega,absorbing_s; version = "BlockAc", blackBox = "mg", intergridType = "BI", relaxType = "Jacobi", levels = 2, recursive_calls = 2, nu = [1;2S])

    direct_setup = [];
    mg_setup = [];

    n = size(mu)
    dim = length(n)

    A_s = blockdiag(AsBlocks...)
    A1_s = AsBlocks[1]; 
    A2_s = AsBlocks[2];
    if dim == 3 
        A3_s = AsBlocks[3]
    end
    Hp,Ap = getPrecBlocks(A_s,Bt,B,C,mu,rho,omega,absorbing_s;version);

    if blackBox == "direct"

        LUA1s = lu(A1_s);
        LUA2s = lu(A2_s);
        if dim == 3
            LUA3s = lu(A3_s)
        end

        if version == "BlockAc"
            LUHp = lu(Hp);
        elseif version == "BFBt" || version == "Fp"
            LUHp = lu(B*Bt); # abuse of notation alert
        end

        direct_setup = [[LUA1s];[LUA2s];[LUHp]];

    elseif blackBox == "mg"

        recursive_calls = 2;
        if dim == 2
            relaxParam = [0.8;0.8;0.3];
            centers = [false;false];
            faces1 = [true;false];
            faces2 = [false;true];
        elseif dim == 3
            relaxParam = [0.8;0.8;0.2];
            centers = [false;false;false];
            faces1 = [true;false;false];
            faces2 = [false;true;false];
            faces3 = [false;false;true];
            nu = [2;2];
            intergridType = "mixed";
            R_A3,P_A3,Ac_A3,LUAcA3 = myMGsetup(A3_s,n,levels,faces3;intergridType);
        end

        R_Hp,P_Hp,Ac_Hp,LUAcHp = myMGsetup(Hp,n,levels,centers;intergridType);
        R_A1,P_A1,Ac_A1,LUAcA1 = myMGsetup(A1_s,n,levels,faces1;intergridType);
        R_A2,P_A2,Ac_A2,LUAcA2 = myMGsetup(A2_s,n,levels,faces2;intergridType);

        Rs = [[R_A1];[R_A2];[R_Hp]];
        Ps = [[P_A1];[P_A2];[P_Hp]];
        Acs = [[Ac_A1];[Ac_A2];[Ac_Hp]];
        LUAcs = [[LUAcA1];[LUAcA2];[LUAcHp]];
        if dim == 3
            Rs = [[R_A1];[R_A2];[R_A3];[R_Hp]];
            Ps = [[P_A1];[P_A2];[P_A3];[P_Hp]];
            Acs = [[Ac_A1];[Ac_A2];[Ac_A3];[Ac_Hp]];
            LUAcs = [[LUAcA1];[LUAcA2];[LUAcA3];[LUAcHp]];
        end

        println("~~~~~~~~~~ MG setup completed ~~~~~~~~~~")
        println("relax param is: ",relaxParam)
        println("pre+post relax num is: ",nu[1]+nu[2])
        println("intergrid: ",intergridType)

        mg_setup = [[intergridType];[relaxType];[relaxParam];[nu];[levels];[recursive_calls];[Rs];[Ps];[Acs];[LUAcs]];

    end

    setup = [[direct_setup];[mg_setup]];

    return setup, Hp, Ap
end

