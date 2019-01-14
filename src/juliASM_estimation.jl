"""
    create_Ux(alpha,beta)

Function that generates a function to compute potential energy for
a region with N CpG sites and with parameters α and β.

# Examples
```julia-repl
julia> Ux_fun = create_Ux(N,alpha,beta)
```
"""
function create_Ux(a::Float64, b::Float64)
    function Ux(x::Array{Int64,1})::Float64
        -a * sum(x) - b * sum(x[2:end] .* x[1:(end-1)])
    end
    return Ux
end # end create_Ux
"""
comp_full_stats(xobs)

Compute vector S of sufficient statistics under model with N CpG sites.
This function is only valid for fully observed data since we are summing
over to obtain overall statistic S=[S1,S2] where
S1 = ∑ T1(x); S2 = ∑ T2(x)
from all complete observations.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs=JuliASM.gen_full_data(10,4);
julia> Svec=JuliASM.comp_full_stats(xobs)
2-element Array{Int64,1}:
-4
-6
```
"""
function comp_full_stats(xobs::Array{Array{Int64,1},1})::Array{Int64,1}
    M = size(xobs)[1]
    s1 = sum(sum(xobs))
    s2 = sum([sum(xobs[i][2:end] .* xobs[i][1:(end-1)]) for i in 1:M])
    return [s1,s2]
end # comp_full_stats
"""
comp_Z(N,alpha,beta)

Compute partition function of a model with N CpG cites and with parameters
α and β. This is computed using the closed-form expression provided in
REFERENCE. The transfer-matrix method is used to obtain a computationally
efficient expression without the need of recursive expressions.

# Examples
```julia-repl
julia>JuliASM.comp_Z(4,1.0,1.0)
1149.715905067999
```

"""
function comp_Z(N::Int64, a::Float64, b::Float64)::Float64
    # Eigenvalues
    aux1 = 0.5 * exp(-a-b)
    aux2 = exp(2*(a+b)) + exp(2*b)
    aux3 = sqrt((exp(2*a)-1)^2 * exp(4*b) + 4 * exp(2*a))
    l1 = aux1 * (aux2 - aux3)
    l2 = aux1 * (aux2 + aux3)
    LN = [l1^(N-1) 0; 0 l2^(N-1)]

    # Matrix S
    aux1 = 0.5 * exp(-a)
    aux2 = exp(2*(a+b)) - exp(2*b)
    s1 = [-aux1 * (aux2 + aux3); 1]
    s2 = [aux1 * (-aux2 + aux3); 1]
    S = [s1 s2]

    # Matrix Sinv
    s1inv = exp(a) / aux3 * [ -1; 1]
    s2inv = [(-aux2+aux3)/ (2*aux3); 0.5 + 0.5 * aux2 / aux3 ]
    Sinv = [s1inv s2inv]

    # Return a Z>0
    u = [exp(-a/2); exp(a/2)]
    return max(1e-100, u' * S * LN * Sinv * u)
end
"""
comp_scal_fac(k,alpha,beta,alphap)

Compute scaling factor in a model with k unobserved CpG cites with parameters
α and β. This is computed using the closed-form expression provided in
REFERENCE. The transfer-matrix method is used to obtain a computationally
efficient expression without the need of recursive expressions.

# Examples
```julia-repl
julia>JuliASM.comp_scal_fac(1,1.0,1.0,1.0)
4.7048192304864935
```

"""
function comp_scal_fac(k::Int64, a::Float64, b::Float64, ap::Float64)::Float64
    # Return one if k=0
    k==0 && return 1

    # Eigenvalues
    aux1 = 0.5 * exp(-a-b)
    aux2 = exp(2*(a+b)) + exp(2*b)
    aux3 = sqrt((exp(2*a)-1)^2 * exp(4*b) + 4 * exp(2*a))
    l1 = aux1 * (aux2 - aux3)
    l2 = aux1 * (aux2 + aux3)
    LN = [l1^(k-1) 0; 0 l2^(k-1)]

    # Matrix S
    aux1 = 0.5 * exp(-a)
    aux2 = exp(2*(a+b)) - exp(2*b)
    s1 = [-aux1 * (aux2 + aux3); 1]
    s2 = [aux1 * (-aux2 + aux3); 1]
    S = [s1 s2]

    # Matrix Sinv
    s1inv = exp(a) / aux3 * [ -1; 1]
    s2inv = [(-aux2+aux3)/ (2*aux3); 0.5 + 0.5 * aux2 / aux3 ]
    Sinv = [s1inv s2inv]

    # Return scaling factor
    vs = [exp(-ap); exp(ap)]
    ve = [exp(-a/2); exp(a/2)]
    return vs' * S * LN * Sinv * ve
end
"""
comp_lkhd(x,alpha,beta)

Compute likelihood of a partial/full observation x that can be missing
values on both sides (given by 0's). Returns -1 if x is not a valid
observation.

# Examples
```julia-repl
julia> JuliASM.comp_lkhd([0;1;-1;1;0],2.0,2.0)
6.144020836828058e-6
```

"""
function comp_lkhd(x::Array{Int64,1}, a::Float64, b::Float64)::Float64
    # Find kl and kr
    ind_dat = findall(!isequal(0), x)
    kl = ind_dat[1]-1
    kr = length(x)-ind_dat[end]

    # Check if observation are contiguous
    if !all(ind_dat[2:end]-ind_dat[1:(length(ind_dat)-1)].==1)
        # return error
        println(stderr,"[$(now())]: Observation $(x) is not contiguous.")
        println(stderr,"[$(now())]: Exiting JuliASM ...")
        exit(1)
    end

    # Get partition function and energy function
    N = length(x)
    Ux = create_Ux(a, b)
    Z = comp_Z(N, a, b)

    # Get scaling factors due to missing data
    sf = comp_scal_fac(kl, a, b, x[ind_dat[1]]*b+a/2) *
         comp_scal_fac(kr, a, b, x[ind_dat[end]]*b+a/2)

    # Return energy function evaluated at x properly scaled
    return exp(-Ux(x[ind_dat])) * sf / Z
end
"""
create_Llkhd(x,alpha,beta)

Create function to compute the minus log-likelihood function for a region with N
CpG sites given the M partial observations.

# Examples
```julia-repl
julia>xobs=JuliASM.gen_full_data(10,4);
julia>LogLike = JuliASM.create_Llkhd(xobs)
```
"""
function create_Llkhd(xobs::Array{Array{Int64,1},1})
    function Llkhd_fun(eta::Array{Float64,1})::Float64
        aux = 0
        a = eta[1]
        b = eta[2]
        for x in xobs
            ind_dat = findall(!isequal(0), x)
            kl = ind_dat[1]-1
            kr = length(x)-ind_dat[end]

            # Check if observation are contiguous
            if !all(ind_dat[2:end]-ind_dat[1:(length(ind_dat)-1)].==1)
                # return error
                println(stderr,"[$(now())]: Observation $(x) is not contiguous.")
                println(stderr,"[$(now())]: Exiting JuliASM ...")
                exit(1)
            end

            # Get partition function and energy function
            N = length(x)
            Ux = create_Ux(a, b)
            Z = comp_Z(N, a, b)

            # Get scaling factors due to missing data
            sf = comp_scal_fac(kl, a, b, x[ind_dat[1]]*b+a/2) *
                 comp_scal_fac(kr, a, b, x[ind_dat[end]]*b+a/2)

            # Add log of it to Llkhd
            aux += -Ux(x[ind_dat]) + log(sf) - log(Z)
        end
        -aux
    end
    return Llkhd_fun
end # end create_Llkhd
"""
est_eta(xobs)

Estimate parameter vector eta=[alpha; beta] based on full or partial
observations. If the number of full observations is smaller than 30,
then a local optimization procedure is used. Otherwise, a global
optimizer is used.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs=JuliASM.gen_full_data(100,4);
julia> JuliASM.est_eta(xobs)
2-element Array{Float64,1}:
-0.021471185237893084
-0.04713465466621405
```

"""
function est_eta(xobs::Array{Array{Int64,1},1})::Array{Float64,1}
    # Local optimizer if fully observed, multiple locals if partially observed
    xfull = filter!(x->!(0 in x),xobs)
    etahat = [NaN, NaN]
    if length(xfull) >= 10
        L = create_Llkhd(xfull)
        eta0 = [0.0, 0.0]
        # NelderMead can't be boxed
        optim = optimize(L, eta0, NelderMead(), Optim.Options(iterations=100);
            autodiff = :forward)
        etahat = optim.minimizer
    else
        # Try different initializations w/ NelderMead
        L = create_Llkhd(xobs)
        aux = Inf
        etas = [[.0,.0],[-.5,-.5],[.5,-.5],[-.5,.5],[.5,.5]]
        for eta0 in etas
            optim = optimize(L, eta0, NelderMead(),
                Optim.Options(iterations=1000); autodiff = :forward)
            if optim.minimum < aux
                aux = optim.minimum
                etahat = optim.minimizer
            end
        end
        # Boxed Simulated Annealing (SAMI)
        # lower = [-2.0, -2.0]
        # upper = [2.0, 2.0]
        # optim = optimize(L, lower, upper, eta0, SAMIN(rt=0.5; f_tol=1e-10,
        #     x_tol=1e-2), Optim.Options(iterations=10^6))
        # etahat = optim.minimizer
        # # BlackBoxOptim: global optimization package.
        # L = create_Llkhd(xobs)
        # optim = bboptimize(L; SearchRange = (-2.0, 2.0), NumDimensions = 2,
        #     MaxTime = 5.0, TraceMode=:silent)
        # # Return estimate vector
        # return best_candidate(optim)
    end
    # Return estimate
    return etahat
end # end est_eta
"""
mle_mult(xobs)

Estimate parameter vector p=[p_1,...,p_{2^N}] based on a multinomial model
and full observations.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs=JuliASM.gen_mult_full_data(100);
julia> JuliASM.mle_mult(xobs)
```

"""
function mle_mult(xobs::Array{Array{Int64,1},1})::Array{Float64,1}

    # Get xcal
    N = size(xobs[1])[1]
    xcal = generate_xcal(N)

    # Get empirical estimates
    phat = zeros(Float64,2^N)
    for x in xobs
        phat[findfirst(y-> x==y, xcal)] += 1
    end
    phat = phat ./ sum(phat)

    # Return phat
    return phat
end # end mle_mult
"""
mle_bin_semi(xobs)

Estimate parameter vector p=[p_1,...,p_{2^N}] based on a binomial model that
assumes a different probability of methylation at each CpG site, as well as
fully observed vectors. This results in a semi-parametric model that grows
in size proportionally to the number of CpG sites.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs=JuliASM.gen_mult_full_data(100);
julia> JuliASM.mle_bin_semi(xobs)
```

"""
function mle_bin_semi(xobs::Array{Array{Int64,1},1})::Array{Float64,1}

    # Get xcal
    N = size(xobs[1])[1]
    M = length(xobs)
    xcal = generate_xcal(N)

    # Get empirical estimate
    θ = zeros(Float64,N)
    [θ[i] = length(findall(x->x[i]==1,xobs)) / M for i in 1:N]

    # Compute probabilities of each pattern
    i = 1
    phat = zeros(Float64,2^N)
    for x in xcal
        phat[i] =prod([θ[j]^isequal(x[j],1) * (1-θ[j])^(1-isequal(x[j],1))
            for j in 1:N])
        i += 1
    end

    # Return phat
    return phat #/ sum(phat)
end # end mle_bin_semi
"""
mle_bin_param(xobs)

Estimate parameter vector p=[p_1,...,p_{2^N}] based on a binomial model
that assumes that each CpG site has the same probability of being methylated,
as well as fully observed vectors. This results in a parametric model that does
not grow with the number of CpG sites considered.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs=JuliASM.gen_mult_full_data(100);
julia> JuliASM.mle_bin_param(xobs)
```

"""
function mle_bin_param(xobs::Array{Array{Int64,1},1})::Array{Float64,1}

    # Get xcal
    N = size(xobs[1])[1]
    M = length(xobs)
    xcal = generate_xcal(N)

    # Get empirical estimates θ=[θ_1,...,θ_N]
    θ = sum([length(findall(x->x==1,x)) / N for x in xobs]) / M

    # Compute probabilities of each pattern
    i = 1
    phat = zeros(Float64,2^N)
    for x in xcal
        k = length(findall(y->y==1,x))
        phat[i] = θ^(k) * (1-θ)^(N-k)
        i+=1
    end

    # Return phat
    return phat #/ sum(phat)
end # end mle_bin_param
