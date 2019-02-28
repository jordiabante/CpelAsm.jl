# Import all package functions
using Test
using JuliASM

# Get test path
JuliASM_path = dirname(pathof(JuliASM))
JuliASM_path = splitdir(JuliASM_path)[1]
test_path = JuliASM_path * "/test/"

# Test set for estimation related functions
@testset "Estimation" begin
    # Check proper computation partition function
    @test comp_Z([4],[1.0],1.0)≈1149.716 atol=1e-3
    # Check proper computation scaling factor
    @test comp_g([1],[1.0],1.0,1.0,1.0)≈3.086 atol=1e-3
    # Check proper likelihood
    @test comp_lkhd([0,1,1,1,0],[5],[1.0],1.0)≈0.986 atol=1e-3
end

# Test set for simulations related functions
@testset "Simulations" begin
    # Check proper xcal
    @test generate_xcal(1)==[[-1],[1]]
end

# Test set for information theory related functions
@testset "Output" begin
    # Check proper E[X] computation
    @test sum(comp_ex([4],[0.0],0.0))≈0.0 atol=1e-3
    # Check proper E[XX] computation
    @test sum(comp_exx([4],[0.0],1.0))≈2.285 atol=1e-3
    # Check proper MML computation
    @test comp_mml(comp_ex([4],[0.0],0.0))≈0.5 atol=1e-3
    # Check proper NME computation
    @test comp_nme(trues(4),[4],[0.0],0.0,zeros(4),zeros(3))≈1.0 atol=1e-3
    # Check proper NMI computation
    @test comp_nmi(1.0,0.0,1.0)≈0.5 atol=1e-3
end
