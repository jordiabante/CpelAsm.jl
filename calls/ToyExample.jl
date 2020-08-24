# Load environment
using Pkg
Pkg.activate(".")

# Load package
using CpelAsm
using Printf

# Define I/O
dir = pwd()
b1 = "$(dir)/test/bam/example.a1.bam"
b2 = "$(dir)/test/bam/example.a2.bam"
fa = "$(dir)/test/fasta/n-masked/example.fa"
vcf = "$(dir)/test/vcf/example.vcf"
out = "$(dir)/test/out/"

# Clean output dir
rm(out,recursive=true)

# Run CpelAsm
CpelAsm.print_log("Running ASM analysis. This should take <2 minutes.")
run_analysis(b1,b2,b1,vcf,fa,out;win_exp=10,n_null=1000,cov_ths=5)

# Print out results
CpelAsm.print_log("Tmml P-values:\n")
println("chr\tstart\tend\tTmml\tP-value\n***\t*****\t***\t****\t*******")
file_cont = split(String(read("$(out)/example_tmml_pvals.bedGraph")),('\t','\n'))[1:(end-1)]
s1 = @sprintf "1\t%s\t%s\t%.4f\t%.4f" file_cont[2] file_cont[3] round(parse(Float64,file_cont[4]),digits=4) round(parse(Float64,file_cont[5]),digits=4)
s2 = @sprintf "1\t%s\t%s\t%.4f\t%.4f\n" file_cont[7] file_cont[8] round(parse(Float64,file_cont[9]),digits=4) round(parse(Float64,file_cont[10]),digits=4)
println(s1)
println(s2)
CpelAsm.print_log("Tnme P-values:\n")
println("chr\tstart\tend\tTnme\tP-value\n***\t*****\t***\t****\t*******")
file_cont = split(String(read("$(out)/example_tnme_pvals.bedGraph")),('\t','\n'))[1:(end-1)]
s1 = @sprintf "1\t%s\t%s\t%.4f\t%.4f" file_cont[2] file_cont[3] round(parse(Float64,file_cont[4]),digits=4) round(parse(Float64,file_cont[5]),digits=4)
s2 = @sprintf "1\t%s\t%s\t%.4f\t%.4f\n" file_cont[7] file_cont[8] round(parse(Float64,file_cont[9]),digits=4) round(parse(Float64,file_cont[10]),digits=4)
println(s1)
println(s2)
CpelAsm.print_log("Tpdm P-values:\n")
println("chr\tstart\tend\tTnme\tP-value\n***\t*****\t***\t****\t*******")
file_cont = split(String(read("$(out)/example_tpdm_pvals.bedGraph")),('\t','\n'))[1:(end-1)]
s1 = @sprintf "1\t%s\t%s\t%.4f\t%.4f" file_cont[2] file_cont[3] round(parse(Float64,file_cont[4]),digits=4) round(parse(Float64,file_cont[5]),digits=4)
s2 = @sprintf "1\t%s\t%s\t%.4f\t%.4f\n" file_cont[7] file_cont[8] round(parse(Float64,file_cont[9]),digits=4) round(parse(Float64,file_cont[10]),digits=4)
println(s1)
println(s2)
