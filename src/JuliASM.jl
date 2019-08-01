__precompile__()

module JuliASM
###################################################################################################
# DEPENDENCIES
###################################################################################################
using Distributed               # For parallelization
using SharedArrays              # For parallelization
using DelimitedFiles            # For delimited files
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
using MultipleTesting           # For multiple hypothesis testing
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

# Data generating functions
export gen_x_mc                     # Generate a methylation vector using MC approach
export gen_ising_full_data          # Generate a set of full observations using slow approach
export gen_ising_part_data          # Generate a set of partial observations using slow approach

# Different computations
export comp_Z                       # Partition function
export est_alpha                    # Estimate α for N=1 case
export est_theta                    # Estimate θ in general case
export comp_lkhd                    # Compute likelihood
export comp_ex                      # Compute E[X] vector
export comp_exx                     # Compute E[XX] vector
export comp_cov                     # Compute Σ matrix
export comp_mml                     # Compute MML in `comp_tobs`
export comp_mml_∇                   # Compute MML in `comp_tnull (uses ∇logZ)
export comp_nme                     # Compute NME in `comp_tobs`
export comp_nme_∇                   # Compute NME in `comp_tnull` (uses ∇logZ)
export comp_nme_mix_mc              # Compute NME in mixture model using MC technique
export comp_nme_mix_exact           # Compute NME in mixture model recursively (N<17)
export comp_uc                      # Compute UC

# General functions
export gen_gffs                     # Function to generate GFF files
export comp_tobs                    # Compute statistis in haplotypes
export comp_tnull                   # Generate null statistis
export comp_pvals                   # Compute p-values
export run_analysis                 # Run it all

end # module
