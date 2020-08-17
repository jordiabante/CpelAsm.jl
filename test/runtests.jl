# Import all package functions
using Test
using Random
using CpelAsm

# Test set for estimation related functions
@testset "Estimation" begin
    # Initialize all variables
    n=[2,2]; a=[-1.0,1.0]; b=1.0; θ=vcat([a,b]...); ∇logZ=CpelAsm.get_grad_logZ(n,θ);
    ex=CpelAsm.comp_ex(n,a,b); exx=CpelAsm.comp_exx(n,a,b);
    Random.seed!(1234); xobs=[CpelAsm.gen_x_mc(n,a,b) for i=1:50];
    θhat=CpelAsm.est_theta_sa(n,xobs);
    # Check proper computation partition function
    @test CpelAsm.comp_Z(n,a,b)≈235.912 atol=1e-3
    # Check proper computation scaling factor
    @test CpelAsm.comp_g(n,a,b,1.0,1.0)≈143.330 atol=1e-3
    # Check proper likelihood
    @test CpelAsm.comp_lkhd([1,1,0,1,1],n,a,b)≈0.232 atol=1e-3
    # Check proper estimation of θ
    @test sum((θhat-θ).^2)<0.25
end

# Test set for information theory related functions
@testset "Output" begin
    # Initialize all variables
    n=[2,2]; a=[-1.0,1.0]; b=1.0; θ=vcat([a,b]...); ∇logZ=CpelAsm.get_grad_logZ(n,θ);
    ex=CpelAsm.comp_ex(n,a,b); exx=CpelAsm.comp_exx(n,a,b);
    # Check proper E[X] computation
    @test sum(ex)≈0.0 atol=1e-3
    # Check proper E[XX] computation
    @test sum(exx)≈1.285 atol=1e-3
    # Check proper MML computation
    @test comp_mml(trues(4),ex)≈0.5 atol=1e-3
    # Check proper MML computation (2)
    @test CpelAsm.comp_mml_∇(n,∇logZ)≈0.5 atol=1e-3
    # Check proper NME computation
    @test CpelAsm.comp_nme(trues(4),n,a,b,ex,exx)≈0.462 atol=1e-3
    # Check proper NME computation (2)
    h = CpelAsm.comp_nme_∇(n,θ,∇logZ)
    @test h≈0.463 atol=1e-3
    # Check proper NME computation (3)
    @test CpelAsm.comp_nme_mix_exact(trues(4),trues(4),n,n,θ,θ)≈0.463 atol=1e-3
    # Check proper computation of UC
    @test CpelAsm.comp_uc(trues(4),trues(4),n,n,θ,θ,h,h)≈0.0 atol=1e-3
end
