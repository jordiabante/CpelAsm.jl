__precompile__()

module JuliASM
##########################################################################
# TODO JORDI
##########################################################################
# ESTIMATION
# 1. Coverage? Check informME and do simulations.
# STATISTICAL ANALYSIS
# 2. Refine permutation test?
##########################################################################
# DEPENDENCIES
##########################################################################
using Random            # RNG
using Dates             # For showing time
using Distributions     # For sampling RVs
using Plots             # For plots
using LaTeXStrings      # For labels in plots
using ForwardDiff       # For auto-differentiation
using Optim             # For local minimizer
using BioAlignments     # For reading-in BAM/SAM files
using GenomicFeatures   # For BED, GFF3, and bigWig files
using GeneticVariation  # For VCF file
using BioSequences      # For FASTA file
using Combinatorics     # For permutation test
using ProgressMeter     # For tracking progress
# using BlackBoxOptim     # Alternative global minimizer
##########################################################################
# INCLUDES
##########################################################################
include("juliASM_plots.jl")
include("juliASM_estimation.jl")
include("juliASM_simulations.jl")
include("juliASM_bioinformatics.jl")
include("juliASM_output.jl")
##########################################################################
# EXPORTS
##########################################################################
export generate_xcal
export gen_ising_full_data
export gen_ising_part_data
export gen_mult_full_data
export gen_mult_part_data
export comp_Z
export comp_scal_fac
export comp_lkhd
export est_alpha
export est_eta
export perm_test
export mle_mult
export mle_bin_semi
export mle_bin_param
export mle_asymptotics
export comp_estimates
export read_bam
export read_vcf
export read_gff_chr
export comp_ex
export comp_exx
export comp_mml
export comp_shanH
export comp_mi
export run_asm_analysis

end # module
