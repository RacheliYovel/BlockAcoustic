using KrylovMethods
using Multigrid


function solveGMRES_blockprec(ElasticOperator,b,version,blackBox,AsBlocks,Bt,B,Ap,Hp,ns,setup,inner,maxit,tol; solveA2, A2);

    direct_setup = setup[1];
    mg_setup = setup[2];

    n = ns[1]
    nf = ns[2]
    nf1 = ns[3]
    nf2 = ns[4]
    nc = ns[end]
    if length(n) == 3
        nf3 = ns[5]
    end
    
    if blackBox == "direct" # 2D only
        LUA1s = direct_setup[1];
        LUA2s = direct_setup[2];
        LUHp = direct_setup[3];
    elseif blackBox == "mg"
        intergridType = mg_setup[1];
        relaxType = mg_setup[2];
        relaxParam = mg_setup[3];
        nu = mg_setup[4];
        levels = mg_setup[5];
        recursive_calls = mg_setup[6];
        Rs = mg_setup[7];
        Ps = mg_setup[8];
        Acs = mg_setup[9];
        LUAcs = mg_setup[10];
    end

    
    function BlockTriangularPrec(r)
    
        ru = r[1:nf]; 
        rp = r[nf+1:nf+nc];

        if version == "BlockAc"
            rp = B*ru - Ap(rp); 

        elseif version == "Fp" || version == "BFBt"
            ru1 = ru[1:nf1];
            ru2 = ru[nf1+1:nf];
            Ainvru1 = LUA1s\ru1;
            Ainvru2 = LUA2s\ru2;
            Ainvru = [Ainvru1 ; Ainvru2];
            rp = rp - B * Ainvru;

        end

        if blackBox == "direct"
            if version == "Fp" || version == "BFBt"
                Linvrp = LUHp\rp; # abuse of notation alert: LUHp = lu(B*Bt)
                BFBtLinvrp = Hp*Linvrp;
                if version == "Fp"
                    ep = BFBtLinvrp; 
                else
                    ep = LUHp\BFBtLinvrp;
                end

            elseif version == "BlockAc"
                ep = LUHp\rp; 
            end
        elseif blackBox == "mg"
            ep = MGcycle(Hp,rp,0.0*rp,relaxParam,nu[1],nu[2],levels,recursive_calls,Rs[end],Ps[end],Acs[end],LUAcs[end]; relaxType);
        end

        RHS = ru - Bt*ep;
        RHS1 = RHS[1:nf1];
        RHS2 = RHS[nf1+1:nf1+nf2];

        if blackBox == "direct"
            eu1 = LUA1s\RHS1;
            eu2 = LUA2s\RHS2;
        elseif blackBox == "mg"
            eu1 = MGcycle(AsBlocks[1],RHS1,0.0*ru[1:nf1],relaxParam,nu[1],nu[2],levels,recursive_calls,Rs[1],Ps[1],Acs[1],LUAcs[1]; relaxType);
            eu2 = MGcycle(AsBlocks[2],RHS2,0.0*ru[nf1+1:nf1+nf2],relaxParam,nu[1],nu[2],levels,recursive_calls,Rs[2],Ps[2],Acs[2],LUAcs[2]; relaxType);
            if length(n) == 3
                RHS3 = RHS[nf1+nf2+1:end];
                eu3 = MGcycle(AsBlocks[3],RHS3,0.0*ru[nf1+nf2+1:end],relaxParam,nu[1],nu[2],levels,recursive_calls,Rs[3],Ps[3],Acs[3],LUAcs[3]; relaxType);
            end
        end

        eu = [eu1 ; eu2];
        if length(n) == 3
            eu = [eu1 ; eu2 ; eu3];
        end

        return [eu;ep]

    end

    function CSLP(r)
        if length(r) == nf1
            e = MGcycle(AsBlocks[1],r,0.0*r,relaxParam,nu[1],nu[2],levels,recursive_calls,Rs[1],Ps[1],Acs[1],LUAcs[1]; relaxType);
        elseif length(r) == nf2
            e = MGcycle(AsBlocks[2],r,0.0*r,relaxParam,nu[1],nu[2],levels,recursive_calls,Rs[2],Ps[2],Acs[2],LUAcs[2]; relaxType);
        elseif length(r) == nc
            e = MGcycle(Hp,r,0.0*r,relaxParam,nu[1],nu[2],levels,recursive_calls,Rs[end],Ps[end],Acs[end],LUAcs[end]; relaxType);
        end
        return e
    end

    e = fgmres(ElasticOperator, b, inner, maxIter = maxit, M = BlockTriangularPrec, out = 2, tol = tol , flexible = true)[1]
    if solveA2 == true
        println("~~~~~~~~~~~ solution of A2 block: ~~~~~~~~~~~")
        e = fgmres(A2, b[nf1+1:nf1+nf2], inner, maxIter = maxit, M = CSLP, out = 2, tol = tol , flexible = true)[1]
    end
    
end