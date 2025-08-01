function MGcycle(A,b,x,w,nu1,nu2,levels,recursive_calls,R_arr,P_arr,Ac_arr,LUAcoarsest; coarseSolve = "exact", relaxType = "Jacobi")

	if levels == 1
        if coarseSolve == "exact"
		    return LUAcoarsest \ b
        elseif coarseSolve == "gmres3"
            println("coarse solve iterations")
            inner = 3;
            e = fgmres(A, b, inner, maxIter = 1, out = 2, tol = 1e-6 , flexible = true)[1]
            println("finished coarse solve iterations")
            return e
        elseif coarseSolve == "jac"
            tol = 1e-6;
            e,_ = dampedJac(A,b,0.2,(0.0 + 0.0*1im)*zeros(size(b)),tol,10);
            return e
        elseif coarseSolve == "GS"
            tol = 1e-6;
            e,_ = GaussSeidel(A,b,(0.0 + 0.0*1im)*zeros(size(b)),tol,5);
            return e
        end
	end
	
	# pre-smoothing
    tol = 1e-5;
    if relaxType == "Jacobi"
        x,_ = dampedJac(A,b,w[1],x,tol,nu1);
    elseif relaxType == "GS"
        x,_ = GaussSeidel(A,b,x,tol,nu1);
    end


    # compute and restrict the residual
    r = b - A * x;

    R = R_arr[1];
    P = P_arr[1];
    Ac = Ac_arr[1]

    R_arr = R_arr[2:end];
    P_arr = P_arr[2:end];
    Ac_arr = Ac_arr[2:end];

    rc = R * r;

    # solve the error-residual equation directly or recursively
    ec = 0.0 .* rc; # initial guess
	if levels == 2
		recursive_calls = 1;
	end
    for j=1:recursive_calls
        ec = MGcycle(Ac,rc,ec,w[2:end],nu1,nu2,levels-1,recursive_calls,R_arr,P_arr,Ac_arr,LUAcoarsest ; coarseSolve);
    end

    e = P * ec;
    x = x + e;

    if relaxType == "Jacobi"
        x,_ = dampedJac(A,b,w[1],x,tol,nu2);
    elseif relaxType == "GS"
        x,_ = GaussSeidel(A,b,x,tol,nu2);
    end

    return x

end


function dampedJac(A,b,w,x,tol,maxit)

    D = diag(A); # vector
    wDinv = w ./ D;
    r = b - A * x;
    iter = 0;

    for i = 1:maxit
        if (norm(r) > tol)
            x .+= wDinv .* r;
            r = b - A * x;
            iter = iter + 1;
        end
    end

    return x, iter
end

function GaussSeidel(A,b,guess,tol,maxit)

    x = copy(guess);
    r = b - A * x;
    iter = 0;
    L = tril(A)

    for k = 1:maxit
        if (norm(r) > tol)
            x = x + L\r
            r = b - A * x;
            iter = iter + 1;
        end
    end

    return x, iter
end


function intergrid1D(n,isnodal; intergridType = "BI") # Neumann, put the number of cells as n, isnodal is a scalar
    
    if isnodal == true

        if intergridType == "BI" || intergridType == "mixed" # [1 2 1] linear interpolation
            P = spdiagm(-1 => 0.5 .* ones(n), 0 => ones(n+1), 1 => 0.5 .* ones(n));
            P = P[:,1:2:end];
        elseif intergridType == "high" # [1 4 6 4 1] cubic interpolation
            P = spdiagm(-2 => 0.125 .* ones(n-1), -1 => 0.5 .* ones(n), 0 => 0.75 .* ones(n+1), 1 => 0.5 .* ones(n), 2 => 0.125 .* ones(n-1));
            P = P[:,1:2:end];
        end

        R = 0.5 .* P';

    else 
        
        P = spdiagm(-1 => 0.25 .* ones(n-1), 0 => 0.75 .* ones(n), 1 => 0.75 .* ones(n-1), 2 => 0.25 .* ones(n-2));
        P = P[:,2:2:end];

        # keep same ratios in boundaries to get better convergence
        # P[:,1] = 2 .* P[:,1] ./ sum(P[:,1]);
        # P[:,end] = 2 .* P[:,end] ./ sum(P[:,end]);

        # injection in boundaries to fit the choice we did in the monolithic
        P[1,1] = 1;
        P[end,end] = 1;

        if intergridType == "mixed" # [1 1] nearest neighbor interpolation
            R = spdiagm(0 => 0.5 .* ones(n), 1 => 0.5 .* ones(n-1));
            R = R[:,2:2:end];
            R = 1.0 .* R';
        else # [1 3 3 1] quadratic interpolation
            R = 0.5 .* P';
        end

    end
    
    return R,P
end


function intergrid(n,isnodal ; intergridType = "BI") # Neumann, put an array of number of cells as n, isnodal is an array

    dim = length(n);

    R,P = intergrid1D(n[1],isnodal[1] ; intergridType);
    R = [R];
    P = [P];
    for i=2:dim
        temp,_ = intergrid1D(n[i],isnodal[i] ; intergridType);
        R = [[temp] ; R];
        _,temp = intergrid1D(n[i],isnodal[i] ; intergridType);
        P = [[temp] ; P];
    end
    if dim > 1
        R = kron(R...);
        P = kron(P...);
    end
    
    return R,P
end

function myMGsetup(A,n,levels,isnodal;intergridType = "BI",coarseSolve = "exact",CGA = "Galerkin")

    n_arr = (0.5 .^ (0:levels-1))' * 1.0 .* n;
    n_arr = Int.(n_arr);

    R1,P1 = intergrid(n,isnodal ; intergridType);
    Ac1 = R1 * A * P1;

    R_arr = [R1];
    P_arr = [P1];
    Ac_arr = [Ac1];
    if CGA == "Galerkin"
        for i=2:levels-1
            R_temp,P_temp = intergrid(n_arr[:,i],isnodal ; intergridType);
            R_arr = [R_arr ; [R_temp]];
            P_arr = [P_arr ; [P_temp]];
        
            Ac_temp = R_temp * Ac_arr[i-1] * P_temp;
            Ac_arr = [Ac_arr ; [Ac_temp]];
        end
    else Ac_arr = CGA
    end
    
    if coarseSolve == "exact"
        LUAcoarsest = lu(Ac_arr[end]);
    else 
        LUAcoarsest = 0.0;
    end

    return R_arr,P_arr,Ac_arr,LUAcoarsest

end