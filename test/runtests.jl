# Import all package functions
using juliASM

# Required here
using Test
using GenomicFeatures

# Get test path
juliASM_path = dirname(pathof(juliASM))
juliASM_path = splitdir(juliASM_path)[1]
test_path = juliASM_path * "/test/"

# Start testing
println("Testing juliASM...")

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
    @test comp_shanH(4,0.0,0.0)≈4.0 atol=1e-3
    # Check proper Mutual Information computation
    @test comp_mi(4,[-10.0,-10.0],[10.0,10.0])≈1.0 atol=1e-3
end

# Test set for bioinformatics related functions
@testset "Bioinformatics" begin
    gff_path = test_path * "example.gff"
    # Check proper gff reading
    @test GFF3.seqid(juliASM.read_gff_chr(gff_path,"chrTest")[1])=="chrTest"
end
