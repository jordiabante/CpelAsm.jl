###################################################################################################
# FUNCTIONS
###################################################################################################
"""
    `generate_xcal(N)`

Generate state space 𝒳 for the methylation vector with N CpG sites.

# Examples
```julia-repl
julia> xcal=CpelAsm.generate_xcal(2)
4-element Array{Vector{Int64},1}:
 [-1, -1]
 [1, -1]
 [-1, 1]
 [1, 1]
```
"""
function generate_xcal(n::Int64)::Array{Vector{Int64},1}

    # Generate iterative object
    xcal = Vector{Int64}[]
    [push!(xcal,2 .* digits(i, base=2, pad=n) - ones(n)) for i in 0:(2^n-1)]

    # Return 𝒳
    return xcal

end # end generate_xcal
"""
    `gen_ising_full_data(M,N;α,β)`

Generate M FULLY observed reads with [N1,...,NK] CpG sites from energy function
U(X) parametrized by [α1,...,αK] and β.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs=CpelAsm.gen_ising_full_data(10,[4])
10-element Array{Array{Int64,1},1}:
 [-1, 1, -1, -1]
 [1, -1, 1, -1]
 [-1, 1, 1, -1]
 [1, 1, 1, -1]
 [1, 1, 1, -1]
 [-1, 1, -1, 1]
 [-1, 1, -1, 1]
 [-1, 1, -1, 1]
 [-1, -1, 1, 1]
 [1, 1, 1, 1]
```
"""
function gen_ising_full_data(m::Int64,n::Vector{Int64};a=zeros(length(n)),b=0.0)::Array{Vector{Int64},1}

    # Generate iterative object
    xcal = generate_xcal(sum(n))

    # Create function for energy computation
    Ux_fun = create_Ux(n,a,b)

    # Compute probability over state space
    p = zeros(2^sum(n))
    p = [exp(-Ux_fun(xcal[k])) for k in 1:length(xcal)]
    p = p ./ sum(p)

    # Sample full observations from model
    indices = rand(Distributions.Multinomial(m, p))
    indices = vcat([repeat([i],indices[i]) for i in 1:length(indices)]...)

    # Return observations
    return xcal[indices]

end # end gen_ising_full_data
"""
    `gen_ising_part_data(M,N;α,β)`

Generate M partially observed reads with [N1,...,NK] CpG sites from energy function U(X) 
parametrized by [α1,...,αK] and β. The values code for:
    1:  Methylated CpG site.
    -1: Unmethylated CpG site.
    0:  Unobserved CpG site.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs=gen_ising_part_data(10,[4])
10-element Array{Vector{Int64},1}:
[-1, -1, -1, -1]
[1, -1, -1, -1]
[-1, 1, -1, -1]
[-1, -1, 1, -1]
[1, -1, 1, -1]
[1, -1, -1, 1]
[0, 1, -1, 1]
[1, -1, 1, 1]
[1, -1, 1, 1]
[-1, 1, 1, 1]
```
"""
function gen_ising_part_data(m::Int64,n::Vector{Int64};a=zeros(length(n)),b=0.0)::Array{Vector{Int64},1}

    # Get full data
    full_data = gen_ising_full_data(m,n;a=a,b=b)

    # Add missing values on sides
    missleft = rand(Distributions.Binomial(sum(n), 0.01),m)
    missright = rand(Distributions.Binomial(sum(n), 0.01),m)
    for i in 1:m
        full_data[i][1:missleft[i]] = fill(0, missleft[i])
        full_data[i][(sum(n)-missright[i]+1):end] = fill(0, missright[i])
    end

    # Return partial observations
    return full_data

end # end gen_ising_part_data
"""
    `gen_mult_full_data(M;N,p)`

Generate M FULLY observed reads from a multinomial model with vector of probabilities 
`p=[p_1,...,p_{2^N}]`.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs=CpelAsm.gen_mult_full_data(2)
2-element Array{Array{Int64,1},1}:
 [1, 1, -1, 1]
 [1, 1, -1, 1]
```
"""
function gen_mult_full_data(m::Int64;n=4,p=1.0/(2.0^n)*ones(2^n))::Array{Vector{Int64},1}

    # Generate iterative object
    xcal = generate_xcal(n)

    # Sample full observations from model
    indices = rand(Distributions.Multinomial(m, p))
    indices = vcat([repeat([i],indices[i]) for i in 1:length(indices)]...)

    # Return observations
    return xcal[indices]

end # end gen_mult_full_data
"""
    `gen_mult_part_data(M;N,p)`

Generate M partially observed reads with N CpG sites from a multinomial model with vector 
probabilities `p=[p_1,...,p_{2^N}]`. The values code for:
    1:  Methylated CpG site.
    -1: Unmethylated CpG site.
    0:  Unobserved CpG site.

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs=CpelAsm.gen_mult_part_data(10)
10-element Array{Array{Int64,1},1}:
 [0, 0, 0, 0]
 [1, -1, 1, -1]
 [0, 1, 1, -1]
 [1, 1, 1, 0]
 [1, 1, 1, 0]
 [0, 1, -1, 1]
 [0, 1, -1, 1]
 [0, 1, -1, 1]
 [0, -1, 1, 1]
 [1, 1, 1, 1]
```
"""
function gen_mult_part_data(m::Int64;n=4,p=1.0/(2.0^n)*ones(2^n))::Array{Vector{Int64},1}

    # Get full data
    full_data = gen_mult_full_data(m;n=n,p=p)

    # Add missing values on sides
    missleft = rand(Distributions.Binomial(n+1, 0.05),m)
    missright = rand(Distributions.Binomial(n+1, 0.05),m)
    for i in 1:m
        full_data[i][1:missleft[i]] = fill(0, missleft[i])
        full_data[i][(n-missright[i]+1):end] = fill(0, missright[i])
    end

    # Return partial observations
    return full_data

end # end gen_mult_part_data