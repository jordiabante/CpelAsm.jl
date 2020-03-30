###################################################################################################
# CONSTANTS
###################################################################################################
const ETA_MAX_ABS=5.0
const CI_THRESH=1.0
###################################################################################################
# FUNCTIONS
###################################################################################################
"""
    `create_Ux([N1,...,NK],[α1,...,αK],β)`

Function that creates a function to compute potential energy for a region with [N1,...,NK],
parameters [α1,...,αK], and correlation β.

# Examples
```julia-repl
julia> Ux_fun = CpelAsm.create_Ux([2,2],[1.0,1.0],1.0);
julia> Ux_fun([1,1,1,1])
-7.0
```
"""
function create_Ux(n::Vector{Int64},a::Vector{Float64},b::Float64)

    # Define function given n, α, and β
    function Ux(x::Vector{Int64})::Float64
        u = a[1] * sum(x[1:n[1]]) + b * sum(x[1:(end-1)] .* x[2:end])
        @inbounds for i in 2:length(n)
            u += a[i] * sum(x[(sum(n[1:(i-1)])+1):(sum(n[1:i]))])
        end
        return -u
    end

    # Return function
    return Ux

end # end create_Ux
"""
    `get_W(N,α,β)`

Function that returns W^{N-1} computed efficiently via sum of rank-1 matrices
assuming α and β.

# Examples
```julia-repl
julia> CpelAsm.get_W(4,0.0,0.0)
2×2 Array{Float64,2}:
 4.0  4.0
 4.0  4.0
```
"""
function get_W(n::Int64, a::Float64, b::Float64)::Array{Float64,2}

    # if n=1 then return identity
    n==1 && return [1.0 0.0;0.0 1.0]

    # Compute relevant quantities
    exp_a = exp(a)
    exp_b = exp(b)
    cosh_a = 0.5*(exp_a+1.0/exp_a)
    sinh_a = exp_a-cosh_a

    # Eigenvalues
    aux1 = exp_b * cosh_a
    aux2 = sqrt(1.0 + exp_b^4*sinh_a^2)
    lambda1N = (aux1-1.0/exp_b*aux2)^(n-1)
    lambda2N = (aux1+1.0/exp_b*aux2)^(n-1)

    # Eigenvectors
    aux1 = -exp_b^2 * sinh_a
    e1 = [aux1-aux2; 1.0]
    e1 /= sqrt(e1'*e1)
    e2 = [aux1+aux2; 1.0]
    e2 /= sqrt(e2'*e2)

    # Return W^{N-1}
    return e1*e1'*lambda1N + e2*e2'*lambda2N

end
"""
    `get_V([α1,α2],β)`

Function that returns V matrix assuming [α1,α2] and β.

# Examples
```julia-repl
julia> CpelAsm.get_V([1.0,1.0],1.0)
2×2 Array{Float64,2}:
 1.0       0.367879
 0.367879  7.38906
```
"""
function get_V(a::Vector{Float64},b::Float64)::Array{Float64,2}

    # Compute auxiliary vars
    exp_b = exp(b)
    exp_a_p = exp(0.5*sum(a))
    exp_a_m = exp(0.5*(a[2]-a[1]))

    # Return V
    return [exp_b/exp_a_p exp_a_m/exp_b;1.0/(exp_b*exp_a_m) exp_b*exp_a_p]

end
"""
    `get_u(α)`

Function that boundary vector for Z computation.

# Examples
```julia-repl
julia> CpelAsm.get_u(0.0)
2-element Vector{Float64}:
 1.0
 1.0
```
"""
function get_u(a::Float64)::Vector{Float64}

    # Compute auxiliary
    exp_aux = exp(a/2.0)

    # Return u
    return [1.0/exp_aux; exp_aux]

end
"""
    `comp_Z([N1,...,NK],[α1,...,αK],β)`

Compute partition function of a model with [N1,...,NK] CpG cites and with
parameters [α1,...,αK] and β.

# Examples
```julia-repl
julia> CpelAsm.comp_Z([1,1,1],[1.0,1.0,1.0],1.0)
155.37102759254836
```
"""
function comp_Z(n::Vector{Int64},a::Vector{Float64},b::Float64)::Float64

    # Boundary conditions.
    y = get_u(a[1])'*get_W(n[1],a[1],b)
    if length(n)>1
        y *= prod([get_V(a[(i-1):i],b)*get_W(n[i],a[i],b) for i in 2:length(n)])
    end
    y *= get_u(a[end])

    # Return Z
    return max(1e-100,y)

end
"""
    `get_grad_logZ([N1,...,NK],θhat)`

Numerically computes the gradient of the log partition function of a model with [N1,...,NK] 
CpG cites and estimated parameter vector θhat. This is equivalent to computing the expected 
value of sufficient statistis (SS).

# Examples
```julia-repl
julia> CpelAsm.get_grad_logZ([1,1,1],[1.0,-1.0,1.0,0.0])
4-element Array{Float64,1}:
  0.7615941559702503
 -0.7615941559702503
  0.7615941559335818
 -1.1600513167692226
```
"""
function get_grad_logZ(n::Vector{Int64},θhat::Vector{Float64})::Vector{Float64}

    # Define function
    function f(θ::Vector{Float64})
        log(comp_Z(n,θ[1:(end-1)],θ[end]))
    end

    # Return ∇A(θ)
    return Calculus.gradient(f,θhat)

end
"""
    `check_boundary(θhat)`

Function that returns a bool indicating whether model with parameter estimate vector θhat is on the
boundary of the parameter space in any of its components.

# Examples
```julia-repl
julia> CpelAsm.check_boundary([1.0,1.0,1.0,1.0])
false
julia> CpelAsm.check_boundary([1.0,5.0,1.0,1.0])
true
```
"""
function check_boundary(θhat::Vector{Float64})::Bool

    # Return true θhat on boundary.
    return any(isapprox.(abs.(θhat),ETA_MAX_ABS;atol=5e-2)) || any(abs.(θhat).>ETA_MAX_ABS)

end
"""
    `comp_g([R1,...,RK],[α1,...,αK],β,αp1,αp2)`

Compute scaling factor in a model with [R1,...,RK] unobserved CpG cites from each block, with
parameters [α1,...,αK] and β. The transfer-matrix method is used to obtain a computationally
efficient expression without the need of recursive expressions. αp1 and αp2 are determined by the
respective kind of boundary.

# Examples
```julia-repl
julia> CpelAsm.comp_g([1],[1.0],1.0,3.0,3.0)
20.135323991555527
```
"""
function comp_g(r::Vector{Int64},a::Vector{Float64},b::Float64,ap1::Float64,ap2::Float64)::Float64

    # Return one if r=0
    r==0 && return 1.0

    # Boundary conditions
    y = get_u(ap1)'*get_W(r[1],a[1],b)
    if length(r)>1
        y *= prod([get_V(a[(i-1):i],b)*get_W(r[i],a[i],b) for i in 2:length(r)])
    end
    y *= get_u(ap2)

    # Return scaling factor
    return max(1e-100,y)

end
"""
    `comp_lkhd(X,[N1,...,NK],[α1,...,αK],β)`

Compute likelihood of a partial/full observation X that can be missing values anywhere (given by
0's).

# Examples
```julia-repl
julia> CpelAsm.comp_lkhd([1,1,0,1,1],[5],[1.0],1.0)
0.953646032691218
```
"""
function comp_lkhd(x::Vector{Int64},n::Vector{Int64},a::Vector{Float64},b::Float64)::Float64

    # Avoid the rest if N=1
    length(x)==1 && return 0.5 * exp(x[1]*a[1]) / cosh(a[1])

    # Get partition function and energy function
    Z = comp_Z(n,a,b)
    Ux = create_Ux(n,a,b)

    # Find changes to/from 0 in x vector
    zerost = findall(isequal(-1),abs.(x[2:end]) - abs.(x[1:(end-1)])) .+ 1
    zeroend = findall(isequal(1),abs.(x[2:end]) - abs.(x[1:(end-1)])) .+ 1

    # Determine whether it starts/finishes with 0
    x[1]==0 && pushfirst!(zerost,1)
    x[end]==0 && push!(zeroend,length(x)+1)

    # Find subregion label for each CpG site
    subid = vcat([i*ones(Int64,n[i]) for i in 1:length(n)]...)

    # Get overall scaling factor as product of individual factors
    sf = 1.0
    @inbounds for i in 1:length(zerost)
        # Find boundary conditions (other X's or boundary of X)
        ap1 = zerost[i]==1 ? a[1] : 2.0*x[zerost[i]-1]*b+a[subid[zerost[i]]]
        ap2 = zeroend[i]==sum(n)+1 ? a[end] : 2.0*x[zeroend[i]]*b+a[subid[zeroend[i]-1]]

        # Figure out b (block IDs of missing) and r (α indices)
        b_id = subid[zerost[i]:(zeroend[i]-1)]
        n_miss = [count(x->x==n,b_id) for n in unique(b_id)]

        # Call scaling factor function
        sf *= comp_g(n_miss,a[unique(b_id)],b,ap1,ap2)

    end

    # Return energy function properly scaled.
    return exp(-Ux(x)) * sf / Z

end
"""
    `create_Llkhd([N1,...,NK],XOBS)`

Create function to compute the minus log-likelihood function for a region with N CpG sites given
the M partial observations XOBS.

# Examples
```julia-repl
julia> n=[2,2];
julia> xobs=CpelAsm.gen_ising_full_data(20,n);
julia> LogLike=CpelAsm.create_Llkhd(n,xobs);
julia> LogLike([0.0,0.0,0.0])
55.45177444479562
```
"""
function create_Llkhd(n::Vector{Int64},xobs::Array{Vector{Int64},1})

    # Define minus log-likelihood function
    function Llkhd_fun(theta::Vector{Float64})::Float64
        # Get parameters
        aux = 0.0
        a = theta[1:(end-1)]
        b = theta[end]

        # Get energy function and partition function
        Ux = create_Ux(n,a,b)
        logZ = log(comp_Z(n,a,b))

        # Initialize variables used in for loop
        ap1 = ap2 = 0.0                 # Boundary values left and right
        b_id = n_miss = []              # Block IDs of missing CpG; number missing per ID
        zerost = zeroend = [] # Indices of zeros
        subid = vcat([i*ones(Int64,n[i]) for i in 1:length(n)]...)    # ID's of CpGs

        # Contribution of each observation
        for x in xobs

            # Find changes to/from 0 in x vector
            zerost = findall(isequal(-1),abs.(x[2:end]) - abs.(x[1:(end-1)])) .+ 1
            zeroend = findall(isequal(1),abs.(x[2:end]) - abs.(x[1:(end-1)])) .+ 1

            # Determine whether it starts/finishes with 0
            x[1]==0 && pushfirst!(zerost,1)
            x[end]==0 && push!(zeroend,length(x)+1)

            # Get scaling factors due to all missing data
            @inbounds for i in 1:length(zerost)
                # Find boundary conditions (other X's or boundary of X)
                ap1 = zerost[i]==1 ? a[1] : 2.0*x[zerost[i]-1]*b+a[subid[zerost[i]]]
                ap2 = zeroend[i]==sum(n)+1 ? a[end] : 2.0*x[zeroend[i]]*b+a[subid[zeroend[i]-1]]

                # Figure out b (block IDs of missing) and r (α indices)
                b_id = subid[zerost[i]:(zeroend[i]-1)]
                n_miss = [count(x->x==n,b_id) for n in unique(b_id)]

                # Call scaling factor function
                aux += log(comp_g(n_miss,a[unique(b_id)],b,ap1,ap2))
            end

            # Add log of it to Llkhd
            aux += -Ux(x)
        end

        # Return MINUS log-likelihood. Add as many logZ as samples we have
        -aux + length(xobs) * logZ
    end

    # Return function
    return Llkhd_fun

end # end create_Llkhd
"""
    `est_alpha(XOBS)`

Estimate parameter α in N=1 case (first entry in returned vector).

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> n=[1];
julia> xobs=CpelAsm.gen_ising_full_data(100,n);
julia> CpelAsm.est_alpha(xobs)
2-element Array{Float64,1}:
 -0.020002667306849575
  0.0
```
"""
function est_alpha(xobs::Array{Vector{Int64},1})::Vector{Float64}
    # Derivation of estimate:
        # log(p)+log(2) = α-log(cosh(α)); log(1-p)+log(2) = -α-log(cosh(α))
        # α = (log(p)-log(1-p))/2

    # Proportion of X=1
    phat = length(findall(x->x==[1],xobs)) / length(xobs)
    a = 0.5*(log(phat)-log(1.0-phat))

    # Return estimate
    return [min(max(-ETA_MAX_ABS,a),ETA_MAX_ABS),0.0]

end # end est_alpha
"""
    `est_theta_sa([N1,...,NK],XOBS)`

Estimate parameter vector [α1,...,αK,β] using simulated annealing.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> n=[4]
julia> xobs=CpelAsm.gen_ising_full_data(100,n);
julia> CpelAsm.est_theta_sa(n,xobs)
2-element Vector{Float64}:
 -0.051768715945913854
 -0.022598858872278995
```
"""
function est_theta_sa(n::Vector{Int64},xobs::Array{Vector{Int64},1})::Vector{Float64}

    # If N=1, then estimate α
    sum(n)==1 && return est_alpha(xobs)

    # Define boundaries and initialization
    L = create_Llkhd(n,xobs)
    init = zeros(Float64,length(n)+1)
    lower = -ETA_MAX_ABS * ones(Float64,length(n)+1)
    upper = ETA_MAX_ABS * ones(Float64,length(n)+1)
    opts = Optim.Options(iterations=10^6,show_trace=false,store_trace=false)

    # Boxed Simulated Annealing (SAMI)
    optim = Optim.optimize(L,lower,upper,init,SAMIN(rt=1e-4;f_tol=1e-3,verbosity=0),opts)

    # Return estimate
    return optim.minimizer

end # end est_theta_sa
###################################################################################################
# COMPETING MODELS
###################################################################################################
"""
    `mle_mult(XOBS)`

Estimate parameter vector p=[p_1,...,p_{2^N}] based on a multinomial model and full observations.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs=CpelAsm.gen_mult_full_data(100);
julia> CpelAsm.mle_mult(xobs)
16-element Vector{Float64}:
 0.07
 0.08
 0.06
 0.06
 0.08
 0.08
 0.04
 0.05
 0.04
 0.07
 0.07
 0.02
 0.04
 0.12
 0.07
 0.05
```
"""
function mle_mult(xobs::Array{Vector{Int64},1})::Vector{Float64}

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
    `mle_bin_semi(XOBS)`

Estimate parameter vector `p=[p_1,...,p_{N}]`, based on a binomial model that assumes a different
probability of methylation at each CpG site, from potentially partial observations. This results
in a semi-parametric model that grows in size proportionally to the number of CpG sites.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs=CpelAsm.gen_mult_full_data(100);
julia> CpelAsm.mle_bin_semi(xobs)
4-element Array{Float64,1}:
 0.49
 0.51
 0.46
 0.44
```
"""
function mle_bin_semi(xobs::Array{Vector{Int64},1})::Vector{Float64}

    # Get xcal
    N = size(xobs[1])[1]

    # Get empirical estimate
    θ = [length(findall(x->x[i]==1,xobs)) / (length(findall(x->x[i]==1,xobs))+
         length(findall(x->x[i]==-1,xobs))) for i in 1:N]

    # Return phat
    return θ

end # end mle_bin_semi
"""
    `mle_bin_param(XOBS)`

Estimate parameter vector `p=[p_1,...,p_{2^N}]` based on a binomial model
that assumes that each CpG site has the same probability of being methylated,
as well as fully observed vectors. This results in a parametric model that does
not grow with the number of CpG sites considered.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs=CpelAsm.gen_mult_full_data(100);
julia> CpelAsm.mle_bin_param(xobs)
16-element Array{Float64,1}:
 0.07596914062500001
 0.068733984375
 0.068733984375
 0.062187890625
 0.068733984375
 0.062187890625
 0.062187890625
 0.056265234374999994
 0.068733984375
 0.062187890625
 0.062187890625
 0.056265234374999994
 0.062187890625
 0.056265234374999994
 0.056265234374999994
 0.05090664062499999
```
"""
function mle_bin_param(xobs::Array{Vector{Int64},1})::Vector{Float64}

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
    return phat

end # end mle_bin_param
###################################################################################################
# NOT CURRENTLY USED (EM Algorithm)
###################################################################################################
"""
    `keep_region(M,[N1,...,NK],θhat)`

Function that returns a bool indicating whether model with [N1,...,NK] CpG cites and with parameter
estimate vector θhat, from M observations, should be kept for downstream analysis. Not currently 
used.

# Examples
```julia-repl
julia> CpelAsm.keep_region(10,[1,1,1],[1.0,1.0,1.0,1.0])
false
```
"""
function keep_region(m::Int64,n::Vector{Int64},θhat::Vector{Float64})::Bool

    # Compute FIM considering special case N=1
    fim = sum(n)==1 ? 1-(tanh(θhat[1]))^2 : get_fim(n,θhat)
    sum(n)==1 && return 4.0/sqrt(fim*m)<=CI_THRESH

    # Return true if all CIs are below threshold width
    return det(fim)!=0 ? all(4*sqrt.(max.(0.0,diag(inv(fim)))/m).<=CI_THRESH) : false

end
"""
    `get_fim([N1,...,NK],θhat)`

Estimate Fisher information matrix of a model with [N1,...,NK] CpG cites and evaluated at estimated
parameter vector θhat.

# Examples
```julia-repl
julia> CpelAsm.get_fim([1,1,1],[1.0,-1.0,1.0,0.0])
4×4 Array{Float64,2}:
  0.419975    -3.02772e-6  -3.02772e-6  -0.319849
 -3.02772e-6   0.419978     0.0          0.639701
 -3.02772e-6   0.0          0.419975    -0.319852
 -0.319849     0.639701    -0.319852     1.8143
```
"""
function get_fim(n::Vector{Int64},θhat::Vector{Float64})::Array{Float64,2}

    # Define function
    function f(θ::Vector{Float64})
        log(comp_Z(n,θ[1:(end-1)],θ[end]))
    end

    # Return I(θ)
    return Calculus.hessian(f,θhat)

end
"""
    `comp_ES(𝑋,[N1,...,N_K],θ)`

Function that computes expected value of the complete vector of sufficient statistics 𝐄[S(X)|𝑋;θ],
where X is the full methylation state, given the observed methylation vector 𝑋 and parameter
vector θ.

# Examples
```julia-repl
julia> CpelAsm.comp_ES([1,-1,0,-1],[2,2],[0.0,0.0,0.0])
3-element Array{Float64,1}:
  0.0
 -1.0
 -1.0
```
"""
function comp_ES(x::Vector{Int64},n::Vector{Int64},θ::Vector{Float64})::Vector{Float64}

    # Compute marginal of X=x
    p_x = comp_lkhd(x,n,θ[1:(end-1)],θ[end])

    # Compute E[W_n|X;θ]
    w = zeros(Int64,length(x))
    E_W_X = zeros(Float64,length(x))
    @inbounds for i=1:length(x)
        w[i] = 1
        E_W_X[i] = x[i]==0 ? (2.0*comp_lkhd(x+w,n,θ[1:(end-1)],θ[end])/p_x)-1.0 : x[i]
        w[i] = 0
    end

    # Compute E[W_{n}W_{n+1}|X;θ]
    w1 = zeros(Int64,length(x))
    w2 = zeros(Int64,length(x))
    E_WW_X = zeros(Float64,length(x)-1)
    @inbounds for i=1:(length(x)-1)
        if (x[i]!=0)&(x[i+1]!=0)
            E_WW_X[i] += x[i]*x[i+1]
            continue
        end
        if (x[i]!=0)&(x[i+1]==0)
            E_WW_X[i] += x[i]*E_W_X[i+1]
            continue
        end
        if (x[i]==0)&(x[i+1]!=0)
            E_WW_X[i] += E_W_X[i]*x[i+1]
            continue
        end
        if (x[i]==0)&(x[i+1]==0)
            w1[i] = -1
            w2[i+1] = -1
            E_WW_X[i] += comp_lkhd(x+w1+w2,n,θ[1:(end-1)],θ[end])
            w1[i] = 1
            w2[i+1] = 1
            E_WW_X[i] += comp_lkhd(x+w1+w2,n,θ[1:(end-1)],θ[end])
            w1[i] = -1
            E_WW_X[i] += -comp_lkhd(x+w1+w2,n,θ[1:(end-1)],θ[end])
            w1[i] = 0
            w2[i+1] = 0
            E_WW_X[i] *= 2.0/p_x
            E_WW_X[i] += -1.0
        end
    end

    # Get subregion ID for each CpG site
    subid = vcat([i*ones(Int64,n[i]) for i in 1:length(n)]...)

    # Initialize with observed part for k=1,...,K
    Sout = [sum(E_W_X[subid.==id]) for id in unique(subid)]

    # Append K+1 component
    push!(Sout,sum(E_WW_X))

    # Return vector E[S(W)|X;θ]
    return Sout

end # end comp_ES
"""
    `comp_ET(XOBS,[N1,...,N_K],θ)`

Function that computes expectation

    𝐄[T(X_1,...,X_M)|𝑋_1,...,𝑋_M;θ],

where X_i is the complete methylation state that corresponds to the i-th observation, 𝑋_i is the
corresponding observed vector, and θ is the parameter vector the expectation is computed wrt.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> n=[2,2];
julia> xobs=CpelAsm.gen_ising_full_data(100,n);
julia> CpelAsm.comp_ET(xobs,n,[0.0,0.0,0.0])
3-element Array{Float64,1}:
   0.0
 -20.0
  -6.0
```
"""
function comp_ET(xobs::Array{Vector{Int64},1},n::Vector{Int64},θ::Vector{Float64})::Vector{Float64}

    # Return vector E[T({W1,...,Wm})|{X1,...,Xm};θ]
    return sum(map(x -> comp_ES(x,n,θ),xobs))

end # end comp_ET
"""
    `em_alg([N1,...,N_K],XOBS)`

Function that runs a single instance of the EM algorithm given observations in XOBS. Returns a pair
`θhat` and `convergence` flag.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> n=[2,2];
julia> xobs=CpelAsm.gen_ising_full_data(100,n);
julia> CpelAsm.em_alg(n,xobs)
([-0.00114078, -0.10281, -0.0235107], true)
```
"""
function em_alg(n::Vector{Int64},xobs::Array{Vector{Int64},1})::Tuple{Vector{Float64},Bool}

    # Check if we have complete data
    complete = findfirst(x->any(x.==0),xobs)==nothing

    # Initialize θhat
    θhat = rand(Uniform(-1.5,1.5),length(n)+1)

    # Iterate until ||θ-θhat||<ϵ
    ϵ = 5e-3
    em_convergence = false
    @inbounds for i=1:10
        # (Re)-compute E[T|X1,...,XM;θ^{(i-1)}]
        ET = comp_ET(xobs,n,θhat)
        # Set up system
        function f!(U,θ)
            ∇logZ = get_grad_logZ(n,θ)
            @inbounds for j=1:(length(n)+1)
                U[j] =  ∇logZ[j] - ET[j] / length(xobs)
            end
        end
        #  Try to solve NL system
	    sol =
        try
	        # Call solver
	        nlsolve(f!,θhat;iterations=20,ftol=1e-3)
        catch x
	        # Report x error if found
	        print_log("Issue with NL solver.")
	        break
        end
        # Leave if no convergence
        converged(sol) || break
        # Check Euclidean distance
        em_convergence = complete || (sqrt(sum((sol.zero-θhat).^2)/sum(θhat.^2))<=ϵ)
        # Update θhat
        θhat = sol.zero
        # Break if converged
        em_convergence && break
    end

    # Return θhat
    return (θhat,em_convergence)

end # end em_alg
"""
    `est_theta_em([N1,...,N_K],XOBS)`

Function that runs multiple instances of the EM algorithm given observations in XOBS. Returns a pair
(θhat,convergence flag). This is done so as to try different initializations for the EM algorithm.
The θhat that maximizes the log-likelihood is returned as long as at least one instance of the EM
algorithm has converged.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> n=[2,2];
julia> xobs=CpelAsm.gen_ising_full_data(100,n);
julia> CpelAsm.est_theta_em(n,xobs)
([-0.00114078, -0.10281, -0.0235107], true)
```
"""
function est_theta_em(n::Vector{Int64},xobs::Array{Vector{Int64},1})::Tuple{Vector{Float64},Bool}

    # If N=1, then estimate α
    sum(n)==1 && return (est_alpha(xobs),true)

    # Check if we have complete data
    complete = findfirst(x->any(x.==0),xobs)==nothing

    # Initialize
    min_LogLike = Inf
    convergence = false
    LogLike = create_Llkhd(n,xobs)
    θhat = zeros(Float64,length(n)+1)

    # Run EM with a number of different initializations
    @inbounds for i=1:15
        # EM Algorithm: obtain candidate and convergence flag
        θcand,converged = em_alg(n,xobs)
        # If did not converge move onto new initialization
        converged || continue
        # If it converged and minimizes minus log-likelihood, then keep estimate
        if LogLike(θcand)<min_LogLike
            θhat = θcand
            min_LogLike = min_LogLike
        end
        # Set convergence true since at least one init worked
        convergence = true
        # If full data then no need to re-run EM
        complete && break
    end

    # Return θhat
    return (θhat,convergence)

end # end est_theta_em