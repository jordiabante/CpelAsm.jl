__precompile__()

module JuliASM
###################################################################################################
# DEPENDENCIES
###################################################################################################
using Distributed               # For parallelization
using Random                    # For random number generation
using StatsBase                 # For statistics
using Calculus                  # For Hessian estimation
using LinearAlgebra             # For Eigendecomposition of Σ
using Dates                     # For showing time
using Distributions             # For sampling RVs
using Plots                     # For plots
using LaTeXStrings              # For labels in plots
using ForwardDiff               # For auto-differentiation
using Optim                     # For local minimizer
using BioAlignments             # For reading-in BAM/SAM files
using GenomicFeatures           # For BED, GFF3, and bigWig files
using GeneticVariation          # For VCF file
using BioSequences              # For FASTA file
using Combinatorics             # For permutation test
using ProgressMeter             # For tracking progress
###################################################################################################
# INCLUDES
###################################################################################################
include("juliASM_plots.jl")
include("juliASM_estimation.jl")
include("juliASM_simulations.jl")
include("juliASM_bioinformatics.jl")
include("juliASM_output.jl")
###################################################################################################
# EXPORTS
###################################################################################################
export generate_xcal
export gen_x_mc
export gen_ising_full_data
export gen_ising_part_data
export gen_mult_full_data
export gen_mult_part_data
export comp_Z
export comp_g
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
export gen_gffs
export read_gff_chr
export comp_ex
export comp_exx
export comp_cov
export comp_mml
export comp_nme
export comp_nme_mix_mc
export comp_nme_mix_exact
export comp_nmi
export comp_uc
export comp_corr
export comp_evec
export comp_tobs
export comp_tnull
export run_analysis

end # module
