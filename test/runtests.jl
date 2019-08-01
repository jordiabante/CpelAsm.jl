# Import all package functions
using Test
using Random
using JuliASM

# Get test path
JuliASM_path = dirname(pathof(JuliASM))
JuliASM_path = splitdir(JuliASM_path)[1]
test_path = JuliASM_path * "/test/"

# Test set for estimation related functions
@testset "Estimation" begin
    # Initialize all variables
    n=[2,2]; a=[-1.0,1.0]; b=1.0; θ=vcat([a,b]...); ∇logZ=JuliASM.get_grad_logZ(n,θ);
    ex=JuliASM.comp_ex(n,a,b); exx=JuliASM.comp_exx(n,a,b);
    Random.seed!(1234); xobs=[JuliASM.gen_x_mc(n,a,b) for i=1:50]
    # Check proper computation partition function
    @test JuliASM.comp_Z(n,a,b)≈235.912 atol=1e-3
    # Check proper computation scaling factor
    @test JuliASM.comp_g(n,a,b,1.0,1.0)≈143.330 atol=1e-3
    # Check proper likelihood
    @test JuliASM.comp_lkhd([1,1,0,1,1],n,a,b)≈0.232 atol=1e-3
    # Check proper estimation of θ
    @test sum((JuliASM.est_theta(n,xobs)-θ).^2)<0.5
end

# Test set for information theory related functions
@testset "Output" begin
    # Initialize all variables
    n=[2,2]; a=[-1.0,1.0]; b=1.0; θ=vcat([a,b]...); ∇logZ=JuliASM.get_grad_logZ(n,θ);
    ex=JuliASM.comp_ex(n,a,b); exx=JuliASM.comp_exx(n,a,b);
    # Check proper E[X] computation
    @test sum(ex)≈0.0 atol=1e-3
    # Check proper E[XX] computation
    @test sum(exx)≈1.285 atol=1e-3
    # Check proper MML computation
    @test JuliASM.comp_mml(trues(4),ex)≈0.5 atol=1e-3
    # Check proper MML computation (2)
    @test JuliASM.comp_mml_∇(n,∇logZ)≈0.5 atol=1e-3
    # Check proper NME computation
    @test JuliASM.comp_nme(trues(4),n,a,b,ex,exx)≈0.462 atol=1e-3
    # Check proper NME computation (2)
    h = JuliASM.comp_nme_∇(n,θ,∇logZ)
    @test h≈0.463 atol=1e-3
    # Check proper NME computation (3)
    @test JuliASM.comp_nme_mix_exact(trues(4),trues(4),n,n,θ,θ)≈0.463 atol=1e-3
    # Check proper computation of UC
    @test JuliASM.comp_uc(trues(4),trues(4),n,n,θ,θ,h,h)≈0.0 atol=1e-3
end
