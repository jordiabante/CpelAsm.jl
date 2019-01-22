# Import all package functions
using JuliASM
using Test

# Get test path
JuliASM_path = dirname(pathof(JuliASM))
JuliASM_path = splitdir(JuliASM_path)[1]
test_path = JuliASM_path * "/test/"

# Start testing
println("Testing JuliASM...")

# Test set for estimation related functions
@testset "Estimation" begin
    # Check proper computation partition function
    @test comp_Z(4,1.0,1.0)≈1149.716 atol=1e-3
    # Check proper computation scaling factor
    @test comp_scal_fac(1,1.0,1.0,1.0)≈4.705 atol=1e-3
    # Check proper likelihood
    @test comp_lkhd([0;1;1;1;0],1.0,1.0)≈0.986 atol=1e-3
end

# Test set for simulations related functions
@testset "Simulations" begin
    # Check proper xcal
    @test generate_xcal(1)==[[-1],[1]]
end

# Test set for information theory related functions
@testset "Output Quantities" begin
    # Check proper MML computation
    @test comp_mml(4,0.0,0.0)≈0.5 atol=1e-3
    # Check proper Shannon's entropy computation
    @test comp_shanH(4,0.0,0.0)≈1.723 atol=1e-3
    # Check proper Mutual Information computation
    @test comp_mi(4,[-10.0,-10.0],[10.0,10.0])≈1.0 atol=1e-3
end
