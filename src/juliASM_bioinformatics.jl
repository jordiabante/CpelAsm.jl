###############################################################################
# CONSTANTS
###############################################################################
const THRESH_MAPQ = -1      # MAPQ threshold
const THRESH_COV = 10       # Coverage threshold
###############################################################################
# STRUCTS
###############################################################################
struct AlignTemp
    strand::String      # '+' or '-'
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
function get_align_strand(paired_end::Bool,flag1::UInt16,flag2::UInt16)::String
    # In single-end mode, OT and CTOT reads will both receive a FLAG of 0 while
    # OB and CTOB both receive a FLAG of 16. In paired end mode:
    #              Read 1       Read 2
    #
    #   OT:         99          147
    #   OB:         83          163
    #   CTOT:      147           99
    #   CTOB:      163           83

    if !paired_end
        if flag1==0
            s="OT"
        elseif flag1==16
            s="OB"
        else
            println(stderr,"[$(now())]: Single end. Was expecting FLAG 0 or and")
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
            println(stderr, "[$(now())]: Paired end. Unexpected flag combination.")
            println(stderr, "[$(now())]: Expected 99/147, 147/99, 83/163, 163/83.")
            println(stderr, "[$(now())]: Exiting JuliASM ...")
            exit(1)
        end
    end

    # Return strand
    return s
end # end get_align_strand
function clean_records(paired_end::Bool, records::Array{BAM.Record,1})::AllAlignTemps
    # Initialize struct
    out = AllAlignTemps(paired_end,[])

    # Consider single vs. paired
    if !paired_end
        # Taake care of flags in single end case
        out.templates = [AlignTemp(get_align_strand(paired_end, BAM.flag(x),
            UInt16(0)), x, BAM.Record()) for x in records]
        return out
    else
        # Retrieve unique template names
        temp_names = unique([BAM.tempname(x) for x in records])
        for temp_name in temp_names
            # Get records with same template name
            temp_recs = filter(x->BAM.tempname(x)==temp_name,records)
            if length(temp_recs)==2
                # Matching pairs: first read is either 99 or 163
                if BAM.position(temp_recs[1]) < BAM.position(temp_recs[2])
                    s = get_align_strand(paired_end, BAM.flag(temp_recs[1]),
                                         BAM.flag(temp_recs[2]))
                    push!(out.templates,AlignTemp(s,temp_recs[1],temp_recs[2]))
                else
                    s = get_align_strand(paired_end, BAM.flag(temp_recs[2]),
                                         BAM.flag(temp_recs[1]))
                    push!(out.templates,AlignTemp(s,temp_recs[2],temp_recs[1]))
                end
            elseif length(temp_recs)==1
                # No matching pair in (expanded) window. This has the problem
                # that we can't determine the strand on which the methylation
                # call is made.
            else
                # Quit since the template name should be unique for fragment
                println(stderr, "[$(now())]: More than 2 records in BAM file")
                println(stderr, "[$(now())]: with same name: $(temp_name).")
                println(stderr, "[$(now())]: Exiting JuliASM ...")
                exit(1)
            end
        end
    end

    # Return struct
    return out
end
"""
    read_bam_coord(BAM_PATH, chr, feat_chrSt, feat_chrEnd, feat_cpg_pos;
                   paired_end=true)

Function that reads in BAM file in BAM_PATH and returns methylation vectors
for those records that overlap with (1-based) genomic coordinates
chr:chrSt-chrEnd at feat_cpg_pos.

# Examples
```julia-repl
julia> read_bam_coord(BAM_PATH,"chr1",30,80,[40,60];paired_end=false)
```
"""
function read_bam_coord(bam_path::String, chr::String, feat_chrSt::Int64,
    feat_chrEnd::Int64, feat_cpg_pos::Array{Int64,1}, chr_size::Int64;
    paired_end::Bool=true)::Array{Array{Int64,1},1}

    # Number of CpG sites is determined by that in the region
    N = length(feat_cpg_pos)

    # Expand window if necessary
    # Here we expand the window in paired end case so that we can find the
    # pairs in order to tell where the methylation call was made.
    if paired_end
        exp_win = [max(1, feat_chrSt-75), min(chr_size, feat_chrEnd+75)]
    else
        exp_win = [feat_chrSt,feat_chrEnd]
    end

    # Get records overlapping window.
    reader = open(BAM.Reader, bam_path, index=bam_path*".bai")
    records_olap = collect(eachoverlap(reader, chr, exp_win[1]:exp_win[2]))

    # Relevant flags in BAM file (both Bismark and Arioc)
    #   XS: meth calls (https://github.com/FelixKrueger/Bismark/tree/master/Docs)
    #   XS: uniqueness (http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
    filter!(x-> (BAM.ismapped(x)) && (haskey(x,"XM")) && (!haskey(x,"XS")) &&
        (BAM.flag(x) in [0,16,83,99,147,163]) &&
        (BAM.mappingquality(x)>THRESH_MAPQ), records_olap)

    # Organize records
    records_org = clean_records(paired_end,records_olap)
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

        # Cross positions of CpG sites if record contains CpG sites
        record_cpg_pos=[x.offset for x in eachmatch(r"[zZ]",meth_call)].-OFFSET
        if length(record_cpg_pos)>0
            olap_cpgs=findall(x-> x in record_cpg_pos, feat_cpg_pos)
            # If overlapping CpG sites then store and add to xobs
            if length(olap_cpgs) > 0
                x = zeros(Int64,N)
                record_cpg_pos = record_cpg_pos[findall(x-> x in feat_cpg_pos,
                    record_cpg_pos)] .+ OFFSET
                x[olap_cpgs] = reduce(replace, ["Z"=>1, "z"=>-1],
                    init=split(meth_call[record_cpg_pos],""))
                push!(xobs,x)
            # elseif (maximum(record_cpg_pos) >= minimum(feat_cpg_pos)) ⊻
            #     (maximum(feat_cpg_pos)>=minimum(record_cpg_pos))
            #     # CpG site positions are not intercalated
            # else
            #     println(stderr,"[$(now())]: Read/s $(BAM.tempname(record.R1))" *
            #       " with flag $(BAM.flag(record.R1)) has unused CpG sites.")
            end
        end
    end # end loop over fragments sequenced

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
        write(writer,pop!(gff_records))
    end
    close(writer)
    # Return empty vector
    return gff_records
end # end write_gff
"""
next_record(rdr, chr_names, recordIn)

Returns next record in stream if in chr_names, and false if the stream reached
its end.

# Examples
```julia-repl
julia> next_record(rdr, chr_names, recordIn)
```
"""
function next_record(reader_vcf::VCF.Reader, chr_names::Array{String,1}, record::VCF.Record)::VCF.Record
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
end
"""
    read_vcf(OUT_GFF_PATH,FASTA_PATH,VCF_PATH,WINDOW_SIZE)

Function that creates a GFF file containing the heterozygous SNPs along with
the positions of the sorrounding CpG sites within a window of WINDOW_SIZE.

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
    chr_size = chr_sizes[findfirst(x->x==curr_chr, chr_names)]

    # Loop over variants
    while !eof(reader_vcf)
        # Read record and check if chromosome is in FASTA file or stream's end
        read!(reader_vcf, record)
        if !(VCF.chrom(record) in chr_names)
            record = next_record(reader_vcf, chr_names, record)
            !VCF.haspos(record) && break
        end
        # Check we have loaded the right chromosome from FASTA
        if !(curr_chr == VCF.chrom(record))
            curr_chr = VCF.chrom(record)
            reader_fasta = open(FASTA.Reader,fasta_path,index=fasta_path*".fai")
            fasta_record = reader_fasta[curr_chr]
            chr_size = chr_sizes[findfirst(x->x==curr_chr, chr_names)]
            close(reader_fasta)
        end
        # Check that variant is heterozygous: 0/1, 1/0, 0/2, 2/0, 1/2, and 2/1.
        gtp = VCF.genotype(record)[1][1]
        if !(split(gtp,r"[/|]")[1] == split(gtp,r"[/|]")[2])
            # Get sequence of relevant DNA
            wSt = VCF.pos(record) - Int64(window_size/2)
            wEnd = VCF.pos(record) + Int64(window_size/2)
            win = [max(1, wSt), min(chr_size, wEnd)]
            wSeq=convert(String,FASTA.sequence(fasta_record,win[1]:win[2]))

            # CG's: obtain the position of CpG sites on the forward strand
            cpg_pos = map(x->getfield(x,:offset),eachmatch(r"CG",wSeq)).+ wSt.-1

            # Store heterozygous SNP & sorrounding CpG sites pos if existing
            if length(cpg_pos)>0
                push!(gff_records,GFF3.Record("$(VCF.chrom(record))\t.\thSNP\t" *
                    "$(VCF.pos(record))\t$(VCF.pos(record))\t.\t.\t.\t" *
                    "Ref=$(VCF.ref(record));Alt=$(VCF.alt(record));GT=$(gtp);" *
                    "N=$(length(cpg_pos));CpGs=$(cpg_pos)"))
            end
        end
        # Check if gff_records object is too big (8 bytes x record) and dump it
        if sizeof(gff_records) >= 80000
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
    write_bG(BEDGRAPH_PATH, BEDGRAPH_RECORDS)

Function that writes records in BEDGRAPH_RECORDS into BEDGRAPH_PATH.

# Examples
```julia-repl
julia> write_bG(BEDGRAPH_PATH,BEDGRAPH_RECORDS)
```
"""
function write_bG(bG_records::Array{Tuple{String,Int64,Int64,Float64},1},
    bG_path::String)::Array{Tuple{String,Int64,Int64,Float64},1}
    open(bG_path, "a") do f
        while length(bG_records)>0
            write(f,join(string.(collect(pop!(bG_records))),"\t"),"\n")
        end
    end
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
    fasta_path::String, window_size::Int64, out_path::String; paired_end::Bool=true)

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

    # Check for existance of at least an output files
    if isfile.(bG_mml1_path,bG_mml2_path,bG_h1_path,bG_h2_path,bG_mi_path)
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

    # Loop over chromosomes
    for chr in chr_names
        # Get windows pertaining to current chromosome
        println(stderr,"[$(now())]: Processing 🧬  $(chr) ...")
        features_chr = read_gff_chr(out_gff_path,chr)
        chr_size = chr_sizes[findfirst(x->x==chr, chr_names)]

        # Loop over windows in chromosome
        for feat in features_chr
            # Get window of interest and CpG sites positions.
            feat_chrSt = GFF3.seqstart(feat) - Int64(window_size/2)
            feat_chrEnd = GFF3.seqend(feat) + Int64(window_size/2)
            feat_atts = GFF3.attributes(feat)
            n = parse(Int64,feat_atts[4][2][1])
            feat_cpg_pos = [parse(Int64,replace(s, r"[\[\] ]" => "")) for s in
                feat_atts[5][2]] # space in regex in on purpose

            # Get vectors from BAM1 overlapping feature
            xobs1 = read_bam_coord(bam1_path, chr, feat_chrSt, feat_chrEnd,
                feat_cpg_pos, chr_size; paired_end=paired_end)

            # Get vectors from BAM2 overlapping feature
            xobs2 = read_bam_coord(bam2_path, chr, feat_chrSt, feat_chrEnd,
                feat_cpg_pos, chr_size; paired_end=paired_end)

            # Estimate each single-allele model
            if length(xobs1) >= THRESH_COV
                eta1 = est_eta(xobs1); # println(eta1)
                mml1 = comp_mml(n,eta1[1],eta1[2])
                h1 = comp_shanH(n,eta1[1],eta1[2])
                push!(bG_mml1_records,(chr,feat_chrSt,feat_chrEnd,mml1))
                push!(bG_h1_records,(chr,feat_chrSt,feat_chrEnd,h1))
            end
            if length(xobs2) >= THRESH_COV
                eta2 = est_eta(xobs2); # println(eta2)
                mml2 = comp_mml(n,eta2[1],eta2[2])
                h2 = comp_shanH(n,eta2[1],eta2[2])
                push!(bG_mml2_records,(chr,feat_chrSt,feat_chrEnd,mml2))
                push!(bG_h2_records,(chr,feat_chrSt,feat_chrEnd,h2))
            end

            # Compute mutual information
            if (length(xobs1)>= THRESH_COV) && (length(xobs2)>= THRESH_COV)
                mi = comp_mi(n,eta1,eta2)
                push!(bG_mi_records,(chr,feat_chrSt,feat_chrEnd,mi))
            end

        end
        # Add values to respective output BigWig files after each chr
        bG_mml1_records = write_bG(bG_mml1_records, bG_mml1_path)
        bG_mml2_records = write_bG(bG_mml2_records, bG_mml2_path)
        bG_h1_records = write_bG(bG_h1_records, bG_h1_path)
        bG_h2_records = write_bG(bG_h2_records, bG_h2_path)
        bG_mi_records = write_bG(bG_mi_records, bG_mi_path)
    end
    println(stderr,"[$(now())]: Done. 😄")
end # end run_asm_analysis
