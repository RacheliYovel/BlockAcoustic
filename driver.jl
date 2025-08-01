include("getModels.jl")
include("getBlockPrec.jl")
include("solveGMRES_blockprec.jl")
include("multigrid.jl")
include("monolithic.jl")


function getAndSolveElasticHelmholtzBlockPrec(n,model,shift,omega; lambda_factor = 1.0, beta = 1.0, version = "BlockAc", blackBox = "mg", intergridType = "BI", relaxType = "Jacobi", nu = [1;2], CGA = "Galerkin", levels = 2, inner = 5, maxit = 200, tol = 1e-6, solveA2 = false)
    ElasticOperator,AsBlocks,Bt,B,C,ns,mu,rho,absorbing_s,b = getElasticHelmholtzBlocks(n,model,shift,omega; lambda_factor, beta)
    if solveA2 == true
        _,ABlocks,_,_,_,_,_,_,_,_ = getElasticHelmholtzBlocks(n,model,0.0,omega; lambda_factor, beta)
        A2 = ABlocks[2]
    else
        A2 = 0.0
    end
    setup, Hp, Ap = BlockPrecSetup(AsBlocks,Bt,B,C,mu,rho,omega,absorbing_s; version,blackBox,levels,intergridType,relaxType,nu)
    solveGMRES_blockprec(ElasticOperator,b,version,blackBox,AsBlocks,Bt,B,Ap,Hp,ns,setup,inner,maxit,tol; solveA2, A2);
end


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                         Table 1                         ") 
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

models = ["const";"linear";"SEAM";"SEG";"Marmousi"];
n_arr = [[[256;128]];[[400;128]];[[128;256]];[[128;256]];[[544;112]]];
omegas = [3;3.5;0.89;2.83;1.8].*pi;
shifts = [0.0;0.1];
omega_factors = [1.0;0.1];
Table1 = zeros(length(omega_factors)*length(shifts),length(models));

for i = 1:length(models)
    for j = 1:length(shifts)
        for k = 1:length(omega_factors)
            n = n_arr[i];
            _,AsBlocks,Bt,B,C,ns,mu,rho,absorbing_s,_ = getElasticHelmholtzBlocks(n, models[i], shifts[j],omegas[i]*omega_factors[k])
            A = blockdiag(AsBlocks...);
            _,Ap = getPrecBlocks(A,Bt,B,C,mu,rho,omegas[i]*omega_factors[k],absorbing_s)
            Xi = B*A - Ap(B);
            AmXi = (1/sum(abs.(A))) * sum(abs.(Xi));
            Table1[k + length(omega_factors)*(j-1),i] = AmXi
        end
    end
end

for l = 1:length(models)
    print(models[l],"      ")
end
println("")
display(Table1)


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                         Table 2                         ") 
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

lambda_factors = [1;10;100;1000];
n_arr = [[[200;64]];[[400;128]];[[800;256]];[[1600;512]]];
omegas = [1.75;3.5;7;14].*pi;

for i=1:length(n_arr)
    println("=================================================================================")
    println("============ Grid size: ",n_arr[i], " ============")
    for j=1:length(lambda_factors)
        println("~~~~~~~~~~~~~ lambda factor: ", lambda_factors[j], " ~~~~~~~~~~~~~")
        getAndSolveElasticHelmholtzBlockPrec(n_arr[i],"linear",0.0,omegas[i]; blackBox = "direct", lambda_factor = lambda_factors[j], inner = 100)
    end
end


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                        Figure 10                        ") 
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

omega_factors = 0.001:0.1:1;
versions = ["BlockAc","Fp","BFBt"];
lambda_factors = [1;1000];
n = [128;64];
omega = 1.5*pi;

for i=1:length(lambda_factors)
    println("lambda factor is ",lambda_factors[i])
    for j=1:length(versions)
        println("version is ",versions[j])
        for k=1:length(omega_factors)
            println("omega factor is ",omega_factors[k])
            getAndSolveElasticHelmholtzBlockPrec(n,"const",0.0,omega*omega_factors[k]; blackBox = "direct", version = versions[j], lambda_factor = lambda_factors[i])
        end
    end
end


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                         Table 5                         ") 
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

n_arr = [[[400;128]];[[800;256]];[[1600;512]]];
omegas = [3.5;7;14].*pi;
levels_arr = [2;3;4];
shifts_blockac = [0.1;0.2;0.4];
shifts_monolithic = [0.1;0.4;0.5];

for i=1:length(n_arr)
    println("=================================================================================")
    println("============ Grid size: ",n_arr[i], " ============")
    for j=1:length(levels_arr)
        println(levels_arr[j],"-level method")
        println("Block-acoustic with shift ", shifts_blockac[j])
        getAndSolveElasticHelmholtzBlockPrec(n_arr[i],"linear",shifts_blockac[j],omegas[i]; levels = levels_arr[j])
        println("Monolithic with shift ", shifts_monolithic[j])
        getAndSolveElasticHelmholtzMonolithic(n_arr[i],"linear",shifts_monolithic[j],omegas[i]; levels = levels_arr[j])
    end
end


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                         Table 6                         ") 
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

n_arr = [[[544;112]];[[1088;224]];[[2176;448]]];
omegas = [1.8;3.6;7.2].*pi;
levels_arr = [2;3;4];
shifts_blockac = [0.2;0.2;0.5];
shifts_monolithic = [0.1;0.2;0.6];

for i=1:length(n_arr)
    println("=================================================================================")
    println("============ Grid size: ",n_arr[i], " ============")
    for j=1:length(levels_arr)
        println(levels_arr[j],"-level method")
        println("Block-acoustic with shift ", shifts_blockac[j])
        getAndSolveElasticHelmholtzBlockPrec(n_arr[i],"Marmousi",shifts_blockac[j],omegas[i]; levels = levels_arr[j])
        println("Monolithic with shift ", shifts_monolithic[j])
        getAndSolveElasticHelmholtzMonolithic(n_arr[i],"Marmousi",shifts_monolithic[j],omegas[i]; levels = levels_arr[j])
    end
end


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                         Table 7                         ") 
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

n = [400;128];
omega = 3.5*pi;
intergrids = ["BI","mixed"]
smoothers = ["Jacobi";"GS"];
shifts = [0.4;0.5];
nu_arr = [[[1;1]];[[1;2]]];

for i=1:length(intergrids)
    println("=================================================================================")
    println("============ Intergrid: ",intergrids[i], " ============")
    for j=1:length(smoothers)
        println("smoother: ",smoothers[j])
        for k=1:length(nu_arr)
            println(nu_arr[k]," relaxation steps")
            for l=1:length(shifts)
                println("shift: ",shifts[l])
                getAndSolveElasticHelmholtzBlockPrec(n,"linear",shifts[l],omega; nu = nu_arr[k], levels = 4, relaxType = smoothers[j], intergridType = intergrids[i])
            end
        end
    end
end


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                         Table 8                         ") 
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

n_arr = [[[128;256]];[[256;512]];[[512;1024]]];
omegas = [2.76;5.52;11.04].*pi;
betas = [1;2/3];
model = "SEG"
shift = 0.1
levels = 3
solveA2 = true

for i=1:length(n_arr)
    println("=================================================================================")
    println("============ Grid size: ",n_arr[i], " ============")
    for j=1:length(betas)
        println("Discretization with beta = ", betas[j])
        getAndSolveElasticHelmholtzBlockPrec(n_arr[i],model,shift,omegas[i]; levels, solveA2, beta = betas[j])
    end
end


println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("                         Table 9                         ") 
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

models = ["Overthrust","OverthrustUmair"]
n_arr = [[[128;128;40]];[[192;192;60]];[[256;256;80]];[[320;320;100]]];
omegas = [1.36;2.02;2.65;3.30].*pi;
levels = [4;3;2];
shifts = [0.4;0.2;0.1];

for i=1:length(models)
    println("=================================================================================")
    println("============ Model: ",models[i], " ============")
    for j=1:length(n_arr)
        println("Size: ",n_arr[j])
        for k=1:length(levels)
            println("levels: ",levels[k],", shift: ",shifts[k])
            getAndSolveElasticHelmholtzBlockPrec(n_arr[j],models[i],shifts[k],omegas[j]; levels = levels[k])
        end
    end
end

println("all done!")