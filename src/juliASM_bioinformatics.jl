###############################################################################
# CONSTANTS
###############################################################################
const THRESH_MAPQ = 20                      # MAPQ threshold
const THRESH_COV = 10                       # Coverage threshold
const FLAGS_ALLOWED = [0,16,83,99,147,163]  # Flags allowed in BAM records
###############################################################################
# STRUCTS
###############################################################################
struct AlignTemp
    strand::String      # Methylation call strand
    R1::BAM.Record      # First record from left to right
    R2::BAM.Record      # Second record from left to right
end
mutable struct AllAlignTemps
    paired::Bool
    templates::Array{AlignTemp,1}
end
###############################################################################
# FUNCTIONS
###############################################################################
"""
    get_align_strand(PAIRED_END,FLAG1,FLAG2)

Function that returns the strand of methylation call based on Bismark's logic.
In single-end mode, OT and CTOT reads will both receive a FLAG of 0 while
OB and CTOB both receive a FLAG of 16. In paired end mode:

             Read 1       Read 2
  OT:         99          147
  OB:         83          163
  CTOT:      147           99
  CTOB:      163           83

# Examples
```julia-repl
julia> JuliASM.get_align_strand(true,UInt16(99),UInt16(147))
"OT"
```
"""
function get_align_strand(pe::Bool,flag1::UInt16,flag2::UInt16)::String

    # Treat SE and PE separately
    if !pe
        if flag1==0
            s="OT"
        elseif flag1==16
            s="OB"
        else
            println(stderr,"[$(now())]: SE: was expecting FLAG 0 or and")
            println(stderr,"[$(now())]: encountered $(flag1) instead.")
            println(stderr,"[$(now())]: Exiting JuliASM ...")
            exit(1)
        end
    else
        if flag1==99 && flag2==147
            s="OT"
        elseif flag1==83 && flag2==163
            s="OB"
        elseif flag1==147 && flag2==99
            s="CTOT"
        elseif flag1==163 && flag2==83
            s="CTOB"
        else
            println(stderr, "[$(now())]: PE: unexpected flag combination.")
            println(stderr, "[$(now())]: Expected 99/147, 147/99, 83/163, 163/83.")
            println(stderr, "[$(now())]: Exiting JuliASM ...")
            exit(1)
        end
    end

    # Return strand
    return s
end # end get_align_strand
"""
    order_bams(PAIRED_END,RECORDS)

Function that returns an AlignTemp object with R1 as the first record in the
forward strand and with the methylation call strand taken from get_align_strand.

# Examples
```julia-repl
julia> JuliASM.order_bams(true,RECORDS)
```
"""
function order_bams(pe::Bool, records::Array{BAM.Record,1})::AlignTemp
    # Check which record appears first in forward strand
    if BAM.position(records[1])<=BAM.position(records[2])
        s = get_align_strand(pe, BAM.flag(records[1]), BAM.flag(records[2]))
        return AlignTemp(s,records[1],records[2])
    else
        s = get_align_strand(pe, BAM.flag(records[2]), BAM.flag(records[1]))
        return AlignTemp(s,records[2],records[1])
    end
end
"""
    clean_records(PAIRED_END,RECORDS)

Function that takes a set of records and returns an AllAlignTemps object that
contains all the properly aligned templates as an array of AlignTemp, which
contains information about the methylation call strand as well as the relevant
BAM records. In the PE case, R1 corresponds to the BAM record that appears
before in the forward strand.

# Examples
```julia-repl
julia> JuliASM.clean_records(true,RECORDS)
```
"""
function clean_records(pe::Bool, records::Array{BAM.Record,1})::AllAlignTemps

    # Initialize struct
    out = AllAlignTemps(pe,[])

    # Consider PE vs. SE
    if pe
        # Handle PE case, first retrieve unique template names
        temp_names = unique([BAM.tempname(x) for x in records])
        for temp_name in temp_names
            # Get records with same template name
            temp_recs = filter(x->BAM.tempname(x)==temp_name,records)

            # There should only be two records with the same template name
            length(temp_recs)==2 && push!(out.templates,order_bams(pe,temp_recs))
            # If length()==1, then no matching pair in (expanded) window. This
            # has the problem that we can't determine the strand on which the
            # methylation call is made. If length()>2, then skip as well since
            # template names should be unique.
        end
    else
        # Handle SE case
        out.templates = [AlignTemp(get_align_strand(pe, BAM.flag(x), UInt16(0)),
                         x, BAM.Record()) for x in records]
    end

    # Return struct
    return out
end # end clean_records
"""
    read_bam_coord(BAM_PATH, CHR, featSt, featEnd, feat_cpg_pos; pe=true)

Function that reads in BAM file in BAM_PATH and returns methylation vectors
for those records that overlap with (1-based) genomic coordinates
chr:chrSt-chrEnd at feat_cpg_pos. The information was taken from:

  XS: meth calls (https://github.com/FelixKrueger/Bismark/tree/master/Docs)
  XS: uniqueness (http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)

# Examples
```julia-repl
julia> read_bam_coord(BAM_PATH,"chr1",30,80,[40,60];pe=false)
```
"""
function read_bam_coord(bam_path::String, chr::String, featSt::Int64,
                        featEnd::Int64, feat_cpg_pos::Array{Int64,1},
                        chr_size::Int64; pe::Bool=true)::Array{Array{Int64,1},1}

    # Number of CpG sites is determined by that in the region
    N = length(feat_cpg_pos)

    # Expand window in PE case to include pair even if it's outside the window
    pe ? exp_win=[max(1,featSt-75),min(chr_size,featEnd+75)] : exp_win=[featSt,featEnd]

    # Get records overlapping window.
    reader = open(BAM.Reader, bam_path, index=bam_path*".bai")
    records_olap = collect(eachoverlap(reader, chr, exp_win[1]:exp_win[2]))

    # Relevant flags in BAM file (both Bismark and Arioc)
    filter!(x-> (BAM.ismapped(x)) && (haskey(x,"XM")) && (!haskey(x,"XS")) &&
           (BAM.flag(x) in FLAGS_ALLOWED) && (BAM.mappingquality(x)>THRESH_MAPQ),
           records_olap)

    # Organize records
    records_org = clean_records(pe,records_olap)
    close(reader)

    # Loop over organized records
    xobs = Array{Int64,1}[]
    for record in records_org.templates
        # Get methylation call and offset (depends on strand where call is made)
        if !(records_org.paired)
            # Obtain methylation call from single end
            meth_call = record.R1[:"XM"]
            record.strand in ["OT","CTOT"] ? OFFSET=1 : OFFSET=2
            OFFSET -= BAM.position(record.R1)
        else
            # Obtain methylation call
            R1_call = record.R1[:"XM"]
            R2_call = record.R2[:"XM"]
            dist_st = abs(BAM.position(record.R2) - BAM.position(record.R1))
            meth_call = R1_call * R2_call[(length(R1_call)+1-dist_st):end]
            record.strand in ["OT","CTOT"] ? OFFSET=1 : OFFSET=2
            OFFSET -= BAM.position(record.R1)
        end

        # Cross positions of CpG sites if template contains CpG sites
        temp_cpg_pos = [x.offset for x in eachmatch(r"[zZ]",meth_call)].-OFFSET
        if length(temp_cpg_pos)>0
            olap_cpgs = findall(x-> x in temp_cpg_pos, feat_cpg_pos)
            # If overlapping CpG sites then store and add to xobs
            if length(olap_cpgs)>0
                x = zeros(Int64,N)
                temp_cpg_pos = temp_cpg_pos[findall(x-> x in feat_cpg_pos,
                                                    temp_cpg_pos)] .+ OFFSET
                x[olap_cpgs] = reduce(replace, ["Z"=>1, "z"=>-1], init=
                                      split(meth_call[temp_cpg_pos],""))
                push!(xobs,x)
            # elseif (maximum(temp_cpg_pos) >= minimum(feat_cpg_pos)) ⊻
            #     (maximum(feat_cpg_pos)>=minimum(temp_cpg_pos))
            #     # CpG site positions are not intercalated
            # else
            #     println(stderr,"[$(now())]: Read/s $(BAM.tempname(record.R1))" *
            #       " with flag $(BAM.flag(record.R1)) has unused CpG sites.")
            end
        end
    end # end loop over templates sequenced

    # Return
    return xobs
end # end read_bam_coord
"""
    write_gff(out_gff_path, gff_records)

Function that appends gff_records into out_gff_path.

# Examples
```julia-repl
julia> write_gff(out_gff_path, gff_records)
```
"""
function write_gff(out_gff_path::String, gff_records::Vector{GFF3.Record})::Vector{GFF3.Record}

    # Append gff records
    output = open(out_gff_path, "a")
    writer = GFF3.Writer(output)
    while length(gff_records)>0
        write(writer,popfirst!(gff_records))
    end
    close(writer)

    # Return empty vector
    return gff_records
end # end write_gff
"""
next_record(READER, CHR_NAMES, RECORD)

Recursive function that returns next record in stream if in chr_names, and false
if the stream reached its end. This function is only used in read_vcf().

# Examples
```julia-repl
julia> next_record(reader, chr_names, record)
```
"""
function next_record(reader_vcf::VCF.Reader, chr_names::Array{String,1},
                     record::VCF.Record)::VCF.Record

    # Check if record in chr_names and end of reader
    if VCF.chrom(record) in chr_names
        # If we're good, return record
        return record
    elseif !eof(reader_vcf)
        # If not end of file, call function again
        read!(reader_vcf, record)
        next_record(reader_vcf, chr_names, record)
    else
        # return empty record if finished stream
        return VCF.Record()
    end
end # end next_record
"""
get_records_ps(READER, RECORD)

Recursive function that returns
This function is only used in read_vcf().

# Examples
```julia-repl
julia> get_records_ps(reader, chr_names, record)
```
"""
function get_records_ps(reader_vcf::VCF.Reader,record::VCF.Record,curr_ps::String,
                        wEnd::Int64)::Tuple{VCF.Record,Int64}

    # Check if record has same PS or end of reader
    if VCF.genotype(record)[1][2]!=curr_ps
        return (record, wEnd)
    elseif !eof(reader_vcf)
        # If not end of file, call function again
        wEnd=VCF.pos(record)
        read!(reader_vcf, record)
        get_records_ps(reader_vcf, record, curr_ps, wEnd)
    else
        # return empty record if finished stream
        return (VCF.Record(), wEnd)
    end

end # end get_records_ps
"""
    read_vcf(OUT_GFF_PATH,FASTA_PATH,VCF_PATH,WINDOW_SIZE)

Function that creates a GFF file containing the heterozygous SNPs along with
the positions of the sorrounding CpG sites within a window of WINDOW_SIZE. The
genotype must be specified using A|B or A/B notations (vertical or diagonal bar).

# Examples
```julia-repl
julia> read_vcf(OUT_GFF_PATH,FASTA_PATH,VCF_PATH,WINDOW_SIZE)
```
"""
function read_vcf(out_gff_path::String, fasta_path::String, vcf_path::String,
                  window_size::Int64)

    # Initialize VCF variables
    gff_records = Vector{GFF3.Record}()
    reader_vcf = open(VCF.Reader, vcf_path)
    record = VCF.Record()

    # Initialize FASTA variables
    reader_fasta = open(FASTA.Reader, fasta_path, index=fasta_path*".fai")
    chr_names = reader_fasta.index.names
    curr_chr = chr_names[1]
    fasta_record = reader_fasta[curr_chr]
    close(reader_fasta)

    # Chrom sizes
    chr_sizes = reader_fasta.index.lengths
    chr_size = chr_sizes[findfirst(x->x==curr_chr,chr_names)]

    # Loop over variants
    while !eof(reader_vcf)

        # Read record if empty and check if chromosome is in FASTA file
        VCF.haspos(record) || read!(reader_vcf, record)
        record = next_record(reader_vcf, chr_names, record)
        VCF.haspos(record) || break

        # Check we have loaded the right chromosome from FASTA
        if !(curr_chr == VCF.chrom(record))
            curr_chr = VCF.chrom(record)
            reader_fasta = open(FASTA.Reader,fasta_path,index=fasta_path*".fai")
            fasta_record = reader_fasta[curr_chr]
            chr_size = chr_sizes[findfirst(x->x==curr_chr, chr_names)]
            close(reader_fasta)
        end

        # Get entire (contiguous) haplotype using the PS field
        wSt = wEnd = VCF.pos(record)
        curr_ps = VCF.genotype(record)[1][2]
        record, wEnd = get_records_ps(reader_vcf, record, curr_ps, wEnd)

        # Obtain DNA sequence from reference and CpG loci from it
        wSt -= Int64(window_size/2)
        wEnd += Int64(window_size/2)
        win = [max(1, wSt), min(chr_size, wEnd)]
        wSeq = convert(String,FASTA.sequence(fasta_record,win[1]:win[2]))
        cpg_pos = map(x->getfield(x,:offset),eachmatch(r"CG",wSeq)).+wSt.-1

        # Store haplotype & corresponding CpG sites if necessary
        length(cpg_pos)>0 &&
        push!(gff_records,GFF3.Record("$(curr_chr)\t.\t$(curr_ps)\t$(win[1])"*
              "\t$(win[2])\t.\t.\t.\tN=$(length(cpg_pos));CpGs=$(cpg_pos)"))

        # Check if gff_records object is too big (8 bytes x record) and dump it
        if sizeof(gff_records)>=80000
            gff_records = write_gff(out_gff_path, gff_records)
        end
    end

    # Dump remaining records
    gff_records = write_gff(out_gff_path, gff_records)

end # end read_vcf
"""
    read_gff_chr(GFF_PATH,chr)

Function that reads in a GFF3 file in GFF_PATH and returns an array of
GenomicFeatures Record objects contained in chromosome chr. The format is the
standard defined by ENSEMBL (https://useast.ensembl.org/info/website/upload/gff3.html).

# Examples
```julia-repl
julia> read_gff_chr(GFF_PATH,"chr1")
```
"""
function read_gff_chr(gff_path::String,chr::String)::Array{GFF3.Record,1}

    # Load genomic features from a GFF3 file
    features = open(collect, GFF3.Reader, gff_path)

    # Keep features in chr
    filter!(x -> GFF3.seqid(x) == chr, features)

    # Return
    return features
end # end read_gff_chr
"""
    plot_gff(GFF_PATH,OUT_PATH)

Function that generates different plots from the JuliASM GFF file.

# Examples
```julia-repl
julia> plot_gff(GFF_PATH,OUT_PATH)
```
"""
function plot_gff(gff_path::String,out_path::String)

    # Initialize GR backend
    Plots.GRBackend()

    # Load genomic features from a GFF3 file
    features = open(collect, GFF3.Reader, gff_path)

    # Plot of number of CpG sites around loci
    x = Array{Int64,1}[]
    x = map(x -> parse(Int64, GFF3.attributes(x)[4][2][1]), features)
    p = histogram(x, yaxis=(:log10,(1,Inf)), bins=:sturges, xticks=0:1:10,
                  fillalpha=0.25, xtickfont=font(12, "Arial"),
                  ytickfont=font(12, "Arial"), ylabel="Counts",
                  xlabel="Number of CpG sites", title=basename(gff_path),
                  legend=:none)
    savefig(p, out_path*basename(gff_path)*".num_cpg_sites.pdf")

    # Plot of distance between loci
    x = []
    prev_st = 1
    curr_chr = unique([GFF3.seqid(x) for x in features])[1]
    for feat in features
        GFF3.seqid(feat)==curr_chr && push!(x,abs(GFF3.seqstart(feat)-prev_st))
        curr_chr = GFF3.seqid(feat)
        prev_st = GFF3.seqstart(feat)
    end
    p = histogram(x, xaxis=(:log10,(1,10^8)),yaxis=(:log10,(1,10^6)), bins=:fd,
                  fillalpha=0.25, xtickfont=font(12, "Arial"),
                  ytickfont=font(12, "Arial"), ylabel="Counts",
                  xlabel="Distance between contiguous CpG sites",
                  title=basename(gff_path), legend=:none)
    savefig(p, out_path*basename(gff_path)*".dist_cpg_sites.pdf")

end # end plot_gff
"""
    write_bG(BEDGRAPH_PATH, BEDGRAPH_RECORDS)

Function that writes records in BEDGRAPH_RECORDS into BEDGRAPH_PATH.

# Examples
```julia-repl
julia> write_bG(BEDGRAPH_PATH,BEDGRAPH_RECORDS)
```
"""
function write_bG(bG_records::Array{Tuple{String,Int64,Int64,Float64},1},
                  bG_path::String)::Array{Tuple{String,Int64,Int64,Float64},1}

    # Open output bedGraph file in append mode
    open(bG_path, "a") do f
        while length(bG_records)>0
            write(f,join(string.(collect(popfirst!(bG_records))),"\t"),"\n")
        end
    end

    # No need to close stream
    return bG_records
end # end write_bG
"""
    run_asm_analysis(BAM1_PATH,BAM2_PATH,VCF_PATH,FASTA_PATH,WINDOW_SIZE,OUT_PATH)

Function that runs juliASM on a pair of BAM files (allele 1, allele 2) given a
VCF file that contains the heterozygous SNPs and a FASTA file that contains the
reference genome.

# Examples
```julia-repl
julia> run_asm_analysis(BAM1_PATH,BAM2_PATH,VCF_PATH,FASTA_PATH,WINDOW_SIZE,OUT_PATH)
```
"""
function run_asm_analysis(bam1_path::String, bam2_path::String, vcf_path::String,
                          fasta_path::String, window_size::Int64, out_path::String;
                          pe::Bool=true)

    # Print initialization of juliASM
    println(stderr,"[$(now())]: Initializing JuliASM ...")

    # Check extension of VCF
    if !isequal(splitext(vcf_path)[2],".vcf")
        println(stderr, "[$(now())]: Wrong extension for VCF input file.")
        println(stderr, "[$(now())]: Exiting JuliASM ...")
        exit(1)
    end

    # Check index file exists as well
    if ! isfile.(bam1_path*".bai", bam2_path*".bai")
        println(stderr, "[$(now())]: At least one BAM index file missing.")
        println(stderr, "[$(now())]: Exiting JuliASM ...")
        exit(1)
    end

    # Check if FASTA index exists
    if ! isfile(fasta_path * ".fai")
        println(stderr, "[$(now())]: FASTA index file not found.")
        println(stderr, "[$(now())]: Exiting JuliASM ...")
        exit(1)
    end

    # Check if output folder exists
    if ! isdir(out_path)
        println(stderr,"[$(now())]: Creating output folder ...")
        mkdir(out_path)
    end

    # BigWig output files
    prefix_sample = split(basename(bam1_path),".")[1]
    bG_mml1_path = "$(out_path)/$(prefix_sample)_mml1.bedGraph"
    bG_mml2_path = "$(out_path)/$(prefix_sample)_mml2.bedGraph"
    bG_h1_path = "$(out_path)/$(prefix_sample)_h1.bedGraph"
    bG_h2_path = "$(out_path)/$(prefix_sample)_h2.bedGraph"
    bG_mi_path = "$(out_path)/$(prefix_sample)_mi.bedGraph"
    bG_pVal_path = "$(out_path)/$(prefix_sample)_pVal.bedGraph"

    # Check for existance of at least an output files
    if isfile.(bG_mml1_path,bG_mml2_path,bG_h1_path,bG_h2_path,bG_mi_path,bG_pVal_path)
        println(stderr, "[$(now())]: At least an output file already exists.")
        println(stderr, "[$(now())]: Make sure you don't overwrite anything.")
        println(stderr, "[$(now())]: Exiting JuliASM ...")
        exit(1)
    end

    # Create gff file with all required information for analysis from vcf
    println(stderr,"[$(now())]: Reading in VCF & FASTA files ...")
    out_gff_path = out_path * split(basename(vcf_path),".")[1] * ".juliasm.gff"
    if isfile(out_gff_path)
        println(stderr, "[$(now())]: Found JuliASM GFF file will be used. If VCF")
        println(stderr, "[$(now())]: file has been updated, change output folder")
        println(stderr, "[$(now())]: or remove existing GFF.")
    else
        println(stderr, "[$(now())]: Generating JuliASM GFF file ... 💪")
        read_vcf(out_gff_path, fasta_path, vcf_path, window_size)
        plot_gff(out_gff_path, out_path)
    end

    # Find chromosomes
    reader_fasta = open(FASTA.Reader, fasta_path, index= fasta_path * ".fai")
    chr_names = reader_fasta.index.names
    chr_sizes = reader_fasta.index.lengths
    chr_list = [(chr_names[i],chr_sizes[i]) for i in 1:length(chr_names)]
    close(reader_fasta)

    # bedGraph records
    bG_mml1_records = Array{Tuple{String,Int64,Int64,Float64},1}()
    bG_mml2_records = Array{Tuple{String,Int64,Int64,Float64},1}()
    bG_h1_records = Array{Tuple{String,Int64,Int64,Float64},1}()
    bG_h2_records = Array{Tuple{String,Int64,Int64,Float64},1}()
    bG_mi_records = Array{Tuple{String,Int64,Int64,Float64},1}()
    bG_pVal_records = Array{Tuple{String,Int64,Int64,Float64},1}()

    # Loop over chromosomes
    n = 0
    tot_feats = 0
    int_feats_1 = 0
    int_feats_2 = 0
    int_feats_mi = 0
    mml1 = mml2 = h1 = h2 = mi = pVal = 0.0
    for chr in chr_names
        # Get windows pertaining to current chromosome
        println(stderr,"[$(now())]: Processing 🧬  $(chr) ...")
        features_chr = read_gff_chr(out_gff_path,chr)
        chr_size = chr_sizes[findfirst(x->x==chr, chr_names)]

        # Loop over windows in chromosome
        for feat in features_chr

            # Get window of interest and CpG sites positions
            tot_feats +=1
            print(stderr,"$tot_feats\r")
            featSt = GFF3.seqstart(feat)
            featEnd = GFF3.seqend(feat)
            feat_atts = GFF3.attributes(feat)
            n = parse(Int64,feat_atts[1][2][1])
            feat_cpg_pos = [parse(Int64,replace(s, r"[\[\] ]" => "")) for s in
                            feat_atts[2][2]] # space in regex in on purpose

            # Get vectors from BAM1 overlapping feature
            xobs1 = read_bam_coord(bam1_path, chr, featSt, featEnd, feat_cpg_pos,
                                   chr_size; pe=pe)
            # println("xobs1:$xobs1")

            # Get vectors from BAM2 overlapping feature
            xobs2 = read_bam_coord(bam2_path, chr, featSt, featEnd, feat_cpg_pos,
                                   chr_size; pe=pe)
            # println("xobs2:$xobs2")

            # Estimate each single-allele model
            if length(xobs1)>=THRESH_COV
                n==1 ? eta1=est_alpha(xobs1) : eta1=est_eta(xobs1)
                # println("eta1:$eta1")
                mml1 = comp_mml(n,eta1[1],eta1[2])
                h1 = comp_shanH(n,eta1[1],eta1[2])
                push!(bG_mml1_records,(chr,featSt,featEnd,mml1))
                push!(bG_h1_records,(chr,featSt,featEnd,h1))
                int_feats_1 += 1
            end
            if length(xobs2)>=THRESH_COV
                n==1 ? eta2=est_alpha(xobs2) : eta2=est_eta(xobs2)
                # println("eta2:$eta2")
                mml2 = comp_mml(n,eta2[1],eta2[2])
                h2 = comp_shanH(n,eta2[1],eta2[2])
                push!(bG_mml2_records,(chr,featSt,featEnd,mml2))
                push!(bG_h2_records,(chr,featSt,featEnd,h2))
                int_feats_2 += 1
            end

            # Compute mutual information
            if (length(xobs1)>=THRESH_COV) && (length(xobs2)>=THRESH_COV)
                mi = comp_mi(n,eta1,eta2)
                push!(bG_mi_records,(chr,featSt,featEnd,mi))
                # Permutation test
                pVal = perm_test(xobs1,xobs2,mi,n)
                push!(bG_pVal_records,(chr,featSt,featEnd,pVal))
                # Increase mi/pVal counter
                int_feats_mi += 1
            end

        end
        # Add values to respective output BigWig files after each chr
        bG_mml1_records = write_bG(bG_mml1_records, bG_mml1_path)
        bG_mml2_records = write_bG(bG_mml2_records, bG_mml2_path)
        bG_h1_records = write_bG(bG_h1_records, bG_h1_path)
        bG_h2_records = write_bG(bG_h2_records, bG_h2_path)
        bG_mi_records = write_bG(bG_mi_records, bG_mi_path)
        bG_pVal_records = write_bG(bG_pVal_records, bG_pVal_path)
    end
    # Print number of features interrogated and finished message
    println(stderr,"[$(now())]: Total number of features: $(tot_feats).")
    println(stderr,"[$(now())]: G1:$(round(100*int_feats_1/tot_feats; sigdigits=2))%.")
    println(stderr,"[$(now())]: G2:$(round(100*int_feats_2/tot_feats;sigdigits=2))%.")
    println(stderr,"[$(now())]: Both:$(round(100*int_feats_mi/tot_feats;sigdigits=2))%.")
    println(stderr,"[$(now())]: Done. 😄")
end # end run_asm_analysis
