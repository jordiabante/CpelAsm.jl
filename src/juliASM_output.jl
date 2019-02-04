"""
`comp_ex(N,α,β)`

Function that computes mean methylation vector E[X] assuming the Ising model
for a methylation state vector of size N and parameters α and β.

# Examples
```julia-repl
julia> JuliASM.comp_ex(4,0.0,0.0)
4-element Array{Float64,1}:
0.5
0.5
0.5
0.5
```
"""
function comp_ex(N::Int64,a::Float64, b::Float64)::Array{Float64,1}

    # Loop over all positions
    x = zeros(Int64,N)
    y = zeros(Float64,N)
    half = rem(N,2)==0 ? div(N,2) : div(N,2)+1
    for j in 1:half
        x[j] = 1
        y[j] = comp_lkhd(x,a,b) - comp_lkhd(-x,a,b)
        x[j] = 0
    end

    # Mirror first half
    y[N:-1:(N-div(N,2)+1)]=y[1:div(N,2)]

    # Return
    return y

end # end comp_ex
"""
`comp_exx(N,α,β)`

Function that computes mean methylation vector E[XiXj] assuming the Ising model
for a methylation state vector of size N and parameters α and β.

# Examples
```julia-repl
julia> JuliASM.comp_exx(4,0.0,0.0)
3-element Array{Float64,1}:
 0.0
 0.0
 0.0
```
"""
function comp_exx(N::Int64,a::Float64, b::Float64)::Array{Float64,1}

    # Loop over all positions
    x = zeros(Int64,N)
    y = zeros(Float64,N-1)
    half = rem(N-1,2)==0 ? div(N-1,2) : div(N-1,2)+1
    for j in 1:half
        x[j] = 1
        x[j+1] = 1
        y[j] = comp_lkhd(x,a,b) + comp_lkhd(-x,a,b)
        x[j] = -1
        x[j+1] = 1
        y[j] -= comp_lkhd(x,a,b) + comp_lkhd(-x,a,b)
        x[j] = 0
        x[j+1] = 0
    end

    # Mirror first half
    y[(N-1):-1:(N-div(N-1,2))]=y[1:div(N-1,2)]

    # Return
    return y

end # end comp_exx
"""
    `comp_mml(ex)`

Function that computes mean methylation level per CpG site assuming the Ising
model for a methylation state vector of size N and parameters α and β.

# Examples
```julia-repl
julia> JuliASM.comp_mml(JuliASM.comp_ex(4,0.0,0.0))
0.5
```
"""
function comp_mml(ex::Array{Float64,1})::Float64

    # Return
    return 0.5/length(ex)*sum(ex) + 0.5

end # end comp_mml
"""
    `comp_shanH(α,β,ex,exx)`

Function that computes normalized Shannon entropy (in nats) assuming the Ising
model for a methylation state vector of size N and parameters α and β.

# Examples
```julia-repl
julia> N=10
julia> a=0.0
julia> b=0.0
julia> JuliASM.comp_shanH(a,b,JuliASM.comp_ex(N,a,b),JuliASM.comp_exx(N,a,b))
2.8906482631788784
```
"""
function comp_shanH(a::Float64,b::Float64,ex::Array{Float64,1},exx::Array{Float64,1})::Float64

    # Return
    return 1.0/log(length(ex)+1)*(log(comp_Z(length(ex),a,b))-a*sum(ex)-b*sum(exx))

end # end comp_shanH
"""
    `comp_mi(h,h1,h2)`

Function that computes the mutual information between haplotype and methylation
state assuming an Ising model for both single-allele and combined models.

# Examples
```julia-repl
julia> JuliASM.comp_mi(1.0,0.5,1.5)
0.0
```
"""
function comp_mi(h::Float64,h1::Float64,h2::Float64)::Float64

    # Return
    return h - 0.5*(h1+h2)

end # end comp_mi
"""
    `bifurcate(xpool,sel)`

Function that divides xpool into two mutually exclusive sets based on sel.

# Examples
```julia-repl
julia> JuliASM.bifurcate([[1,1],[1,-1],[-1,1]], [1,2])
(Array{Int64,1}[[1, 1], [1, -1]], Array{Int64,1}[[-1, 1]])
```
"""
function bifurcate(xpool::Array{Array{Int64,1},1}, sel::Vector{T}) where T <: Integer
    x = xpool[sel]
    asel = trues(length(xpool))
    asel[sel] .= false
    y = xpool[asel]
    return x,y
end
"""
    `perm_test(XOBS1,XOBS2,TOBS,POS_CG)`

# Examples
```julia-repl
julia> Random.seed!(1234);
julia> xobs1=JuliASM.gen_ising_full_data(100,4;a=1.0);
julia> xobs2=JuliASM.gen_ising_full_data(100,4;a=-1.0);
julia> eta1=JuliASM.est_eta(xobs1);
julia> eta2=JuliASM.est_eta(xobs2);
julia> hom=[1,5,10,15];
julia> cpg_pos=[hom,hom,hom];
julia> tobs=JuliASM.comp_mi(cpg_pos,eta1,eta2);
julia> JuliASM.perm_test(xobs1,xobs2,tobs,cpg_pos)
0.0
```
"""
function perm_test(xobs1::Array{Array{Int64,1},1},xobs2::Array{Array{Int64,1},1},
                   tobs::Float64,cpg_pos::Array{Array{Int64,1},1})::Float64
    # Initialize
    better = worse = 0.0
    xpool = vcat(xobs1, xobs2)

    # Loop over combinations and compute statistic
    i = 1
    for subset in combinations(1:length(xpool), length(xobs2))
      # Bifurcate
      test, control = bifurcate(xpool, subset)

      # Estimate parameters
      length(cpg_pos[2])==1 ? eta1=est_alpha(control) : eta1=est_eta(control)
      length(cpg_pos[3])==1 ? eta2=est_alpha(test) : eta2=est_eta(test)

      # Compute MI for partition and compared to that observed
      comp_mi(cpg_pos,eta1,eta2)>tobs ? better+=1.0 : worse+=1.0

      # If enough permutations leave
      i<100 ? i+=1 : break
    end

    # Return p-value
    return better/(worse+better)
end # end perm_test
