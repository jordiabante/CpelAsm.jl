var documenterSearchIndex = {"docs":
[{"location":"main_commands/#CpelAsm","page":"Main Commands","title":"CpelAsm","text":"","category":"section"},{"location":"main_commands/","page":"Main Commands","title":"Main Commands","text":"(Image: CI) (Image: Docs) (Image: MIT license)","category":"page"},{"location":"main_commands/#Generate-GFF-files","page":"Main Commands","title":"Generate GFF files","text":"","category":"section"},{"location":"main_commands/","page":"Main Commands","title":"Main Commands","text":"The first step consists in generating the 2 necessary GFF files: heterozygous and homozygous. The former is the one that contains the haplotypes analyzed, while the latter contains the homozygous regions of the genome.","category":"page"},{"location":"main_commands/","page":"Main Commands","title":"Main Commands","text":"In the following example, the maximum number of CpG sites allowed per haplotype is 25 (n_max=25). In addition, CpelAsm extends the window delimited by the first and last SNP in the haplotype by 100 bp left and right (win_exp=100).","category":"page"},{"location":"main_commands/","page":"Main Commands","title":"Main Commands","text":"# Parameters\nn_max = 25\nwin_exp = 100\n\n# Paths\nproject_path = \"/path/to/SampleName\"\nvcf = \"$(project_path)/vcf/SampleName.phased.vcf\"\nfasta = \"$(project_path)/fasta/SampleName.masked_hg19.fa\"\nhet_gff = \"$(project_path)/CpelAsmOut/SampleName_het.cpelasm.gff\"\nhom_gff = \"$(project_path)/CpelAsmOut/SampleName_hom.cpelasm.gff\"\n\n# Call\ngen_gffs([het_gff,hom_gff],fasta,vcf,win_exp,n_max)\n","category":"page"},{"location":"main_commands/#BedGraphs-MML,-NME,-and-UC","page":"Main Commands","title":"BedGraphs MML, NME, & UC","text":"","category":"section"},{"location":"main_commands/","page":"Main Commands","title":"Main Commands","text":"The next step consists in estimating the allele-specific Ising models in the haplotypes and generating bedGraph files with MML1/2, NME1/2, and UC. CpelAsm can be parallelized when multiple CPUs are available by loading first the package Distributed and then loading CpelAsm through the macro `@everywhere.","category":"page"},{"location":"main_commands/","page":"Main Commands","title":"Main Commands","text":"In the following example, the maximum size of a subregion is fixed to 500 (g_max=500), the minimum average depth is set to 8 (cov_ths=8), and the WGBS reads are trimmed 5 bp on each end (trim=(5,5,5,5)).","category":"page"},{"location":"main_commands/","page":"Main Commands","title":"Main Commands","text":"# Deps\nusing Distributed\n@everywhere using CpelAsm\n\n# Parameters\ng_max = 500\ncov_ths = 8\ntrim = (5,5,5,5)\n\n# Paths\nproject_path = \"/path/to/SampleName\"\nfasta = \"$(project_path)/fasta/SampleName.masked_hg19.fa\"\nhet_gff = \"$(project_path)/cpelasm/SampleName_het.cpelasm.gff\"\nbam1 = \"$(project_path)/bam/SampleName.sort.genome1.bam\"\nbam2 = \"$(project_path)/bam/SampleName.sort.genome2.bam\"\nuc_path = \"$(project_path)/cpelasm/SampleName_uc.bedGraph\"\nmml1_path = \"$(project_path)/cpelasm/SampleName_mml1.bedGraph\"\nmml2_path = \"$(project_path)/cpelasm/SampleName_mml2.bedGraph\"\nnme1_path = \"$(project_path)/cpelasm/SampleName_nme1.bedGraph\"\nnme2_path = \"$(project_path)/cpelasm/SampleName_nme2.bedGraph\"\ntobs_path = [mml1_path,mml2_path,nme1_path,nme2_path,uc_path]\n\n# Call\ncomp_tobs(bam1,bam2,het_gff,fasta,tobs_path;g_max=g_max,cov_ths=cov_ths,trim=trim)\n","category":"page"},{"location":"main_commands/#Generate-Null-Statistics","page":"Main Commands","title":"Generate Null Statistics","text":"","category":"section"},{"location":"main_commands/","page":"Main Commands","title":"Main Commands","text":"The next step consists in generating null statistics to be able to perform hypothesis testing. As discussed in the previous point, CpelAsm can be parallelized when multiple CPUs are available by first loading the package Distributed and then loading CpelAsm through the macro `@everywhere.","category":"page"},{"location":"main_commands/","page":"Main Commands","title":"Main Commands","text":"In the following example, the maximum size of a subregion is fixed to 500 (g_max=500), the minimum average depth is set to 8 (cov_ths=8), and the WGBS reads are trimmed 5 bp on each end (trim=(5,5,5,5)).","category":"page"},{"location":"main_commands/","page":"Main Commands","title":"Main Commands","text":"# Deps\nusing Distributed\n@everywhere using CpelAsm\n\n# Parameters\ncov_ths = 8\ng_max = 500\ntrim = (5,5,5,5)\n\n# Paths\nproject_path = \"/path/to/SampleName\"\nfasta = \"$(project_path)/fasta/SampleName.masked_hg19.fa\"\nbam1 = \"$(project_path)/bam/SampleName.sort.genome1.bam\"\nbam2 = \"$(project_path)/bam/SampleName.sort.genome2.bam\"\nhet_gff = \"$(project_path)/cpelasm/SampleName_het.cpelasm.gff\"\nhom_gff = \"$(project_path)/cpelasm/SampleName_hom.cpelasm.gff\"\nnull_tpdm_path = \"$(project_path)/cpelasm/SampleName_tpdm_null.bedGraph\"\nnull_tmml_path = \"$(project_path)/cpelasm/SampleName_tmml_null.bedGraph\"\nnull_tnme_path = \"$(project_path)/cpelasm/SampleName_tnme_null.bedGraph\"\ntnull_path = [null_tmml_path,null_tnme_path,null_tpdm_path]\nuc_path = \"$(project_path)/cpelasm/SampleName_uc.bedGraph\"\nmml1_path = \"$(project_path)/cpelasm/SampleName_mml1.bedGraph\"\nmml2_path = \"$(project_path)/cpelasm/SampleName_mml2.bedGraph\"\nnme1_path = \"$(project_path)/cpelasm/SampleName_nme1.bedGraph\"\nnme2_path = \"$(project_path)/cpelasm/SampleName_nme2.bedGraph\"\ntobs_path = [mml1_path,mml2_path,nme1_path,nme2_path,uc_path]\n\n# Call\ncomp_tnull(bam,het_gff,hom_gff,fasta,tobs_path,tnull_path;\n    g_max=g_max,cov_ths=cov_ths,trim=trim,n_max=n_max)\n","category":"page"},{"location":"main_commands/#Perform-Allele-Specific-Methylation-Detection","page":"Main Commands","title":"Perform Allele-Specific Methylation Detection","text":"","category":"section"},{"location":"main_commands/","page":"Main Commands","title":"Main Commands","text":"The final step is to compute a p-value for each haplotype for each one of the three quantities (Tmml, Tnme, Tpdm). The returned p-values are corrected using the Benjamini-Hochberg and allow for control of the false discovery rate (FDR).","category":"page"},{"location":"main_commands/","page":"Main Commands","title":"Main Commands","text":"In the following example, hypothesis testing is performed for haplotypes with at most 25 CpG sites (n_max=25), and a minimum of 1,000 null statistics to perform hypothesis testing.","category":"page"},{"location":"main_commands/","page":"Main Commands","title":"Main Commands","text":"# Deps\nusing CpelAsm\n\n# Parameters\nn_max = 25\nn_null = 1000\n\n# Paths\nproject_path = \"/path/to/SampleName\"\nfasta = \"$(project_path)/fasta/SampleName.masked_hg19.fa\"\nbam1 = \"$(project_path)/bam/SampleName.sort.genome1.bam\"\nbam2 = \"$(project_path)/bam/SampleName.sort.genome2.bam\"\nnull_tpdm_path = \"$(project_path)/cpelasm/SampleName_tpdm_null.bedGraph\"\nnull_tmml_path = \"$(project_path)/cpelasm/SampleName_tmml_null.bedGraph\"\nnull_tnme_path = \"$(project_path)/cpelasm/SampleName_tnme_null.bedGraph\"\ntnull_path = [null_tmml_path,null_tnme_path,null_tuc_path]\nuc_path = \"$(project_path)/cpelasm/SampleName_uc.bedGraph\"\nmml1_path = \"$(project_path)/cpelasm/SampleName_mml1.bedGraph\"\nmml2_path = \"$(project_path)/cpelasm/SampleName_mml2.bedGraph\"\nnme1_path = \"$(project_path)/cpelasm/SampleName_nme1.bedGraph\"\nnme2_path = \"$(project_path)/cpelasm/SampleName_nme2.bedGraph\"\ntobs_path = [mml1_path,mml2_path,nme1_path,nme2_path,uc_path]\npVal_tmml_path = \"$(project_path)/cpelasm/SampleName_tmml_pvals.bedGraph\"\npVal_tnme_path = \"$(project_path)/cpelasm/SampleName_tnme_pvals.bedGraph\"\npVal_tpdm_path = \"$(project_path)/cpelasm/SampleName_tpdm_pvals.bedGraph\"\npVal_path = [pVal_tmml_path,pVal_tnme_path,pVal_tuc_path]\n\n# Call\ncomp_pvals(tobs_path,tnull_path,p_path,n_max,n_null)\n","category":"page"},{"location":"main_commands/#Running-CpelAsm-with-a-single-command","page":"Main Commands","title":"Running CpelAsm with a single command","text":"","category":"section"},{"location":"main_commands/","page":"Main Commands","title":"Main Commands","text":"The following example shows how to perform all the steps shown above in a single command.","category":"page"},{"location":"main_commands/","page":"Main Commands","title":"Main Commands","text":"# Deps\nusing Distributed\n@everywhere using CpelAsm\n\n# Parameters\ncov_ths = 8\ng_max = 500\nwin_exp = 100\ntrim = (5,5,5,5)\n\n# Paths\nproject_path = \"/path/to/SampleName\"\nvcf = \"$(project_path)/vcf/SampleName.phased.vcf\"\nfasta = \"$(project_path)/fasta/SampleName.masked_hg19.fa\"\nbam1 = \"$(project_path)/bam/SampleName.sort.genome1.bam\"\nbam2 = \"$(project_path)/bam/SampleName.sort.genome2.bam\"\nbamu = \"$(project_path)/bam/SampleName.sort.unassigned.bam\"\noutdir = \"$(project_path)/cpelasm/\"\n\n# Call\nrun_analysis(bam1,bam2,bamu,vcf,fasta,outdir;g_max=g_max,cov_ths=cov_ths,\n    win_exp=win_exp,trim=trim,n_null=n_null,n_max=n_max)\n","category":"page"},{"location":"toy_example/#CpelAsm","page":"Toy Example","title":"CpelAsm","text":"","category":"section"},{"location":"toy_example/","page":"Toy Example","title":"Toy Example","text":"(Image: CI) (Image: Docs) (Image: MIT license)","category":"page"},{"location":"toy_example/#Toy-Example","page":"Toy Example","title":"Toy Example","text":"","category":"section"},{"location":"toy_example/","page":"Toy Example","title":"Toy Example","text":"The package includes a small toy example, which should take no more than 20 seconds to run, for illustrative purposes. The example in particular consists of 2 haplotypes in a single chromosome of an artificial reference genome. Each haplotype has two alleles: a1 and a2. For simplicity, CpG sites in both alleles are distributed independently according to a Bernoulli distribution. In a1, the probability of methylation is p1=0.8, while in a2 the probability of methylation is p2=0.2. Thus, CpG sites in the former have a mean methylation level (MML) of 0.8, while that of the later is 0.2. Given the symmetry of the problem, however, both alleles have the same Shannon entropy (see Shannon entropy of a flip coin). Thus differential analysis only identifies ASM imbalances in for statistics Tmml and Tpdm. That is, CpelAsm only identifies mean methylation level differences and probability distribution of methylation differences. The output bedGraph files with the suffix pvals contain the results of the statistical test performed for each haplotype. These files can be found in the out/ directory. To run the toy example run the following commands in a julia's REPL:","category":"page"},{"location":"toy_example/","page":"Toy Example","title":"Toy Example","text":"# Load CpelAsm\nusing CpelAsm\n\n# Define I/O\ndir = \"/path/to/CpelAsm.jl/test/\"\nb1 = \"$(dir)/bam/example.a1.bam\"\nb2 = \"$(dir)/bam/example.a2.bam\"\nfa = \"$(dir)/fasta/n-masked/example.fa\"\nvcf = \"$(dir)/vcf/example.vcf\"\nout = \"$(dir)/out/\"\n\n# Run CpelAsm\nrun_analysis(b1,b2,b1,vcf,fa,out;g_max=50,win_exp=10,n_null=50,n_max=10)","category":"page"},{"location":"#CpelAsm","page":"Home","title":"CpelAsm","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: CI) (Image: Docs) (Image: MIT license)","category":"page"},{"location":"#Description","page":"Home","title":"Description","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CpelAsm is a julia package especifically desgined for haplotype allele-specific methylation based on the method in [1]. CpelAsm draws ideas from statistical physics and information theory to detect allele-specific methylation imbalances at the haplotype level.","category":"page"},{"location":"#Testing","page":"Home","title":"Testing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CpelAsm is tested against Julia 1.3.0 on the latest versions of Linux, macOS and Windows.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"#Prerequisites","page":"Home","title":"Prerequisites","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"julia v1.3.0\ngit.","category":"page"},{"location":"#Command","page":"Home","title":"Command","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CpelAsm and dependencies can be installed using the following command:","category":"page"},{"location":"","page":"Home","title":"Home","text":"(v1.3) pkg> add https://github.com/jordiabante/CpelAsm.jl.git","category":"page"},{"location":"#Installation-tests","page":"Home","title":"Installation tests","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"In a julia session run","category":"page"},{"location":"","page":"Home","title":"Home","text":"(v1.3) pkg> test CpelAsm","category":"page"},{"location":"#Authors","page":"Home","title":"Authors","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Jordi Abante","category":"page"},{"location":"#References","page":"Home","title":"References","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"[1] Abante, J., Fang, Y., Feinberg, A.P., Goutsias, J., Detection of haplotype-dependent allele-speciﬁc DNA methylation in WGBS data, Nature Communications 2020 XYZ.","category":"page"}]
}
