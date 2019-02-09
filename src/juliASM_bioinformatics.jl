###############################################################################
# CONSTANTS
###############################################################################
const THRESH_MAPQ = 20                      # MAPQ threshold
const THRESH_COV = 15                       # Coverage threshold
const FLAGS_ALLOWED = [0,16,83,99,147,163]  # Flags allowed in BAM recs
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
    `get_align_strand(PAIRED_END,FLAG1,FLAG2)`

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
    `order_bams(PAIRED_END,RECORDS)`

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
    `clean_records(PAIRED_END,RECORDS)`

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
    `try_olaps(READER,CHR,WINDOW)`

Function that tries to find overlaps in BAM reader @ WINDOW.

# Examples
```julia-repl
julia> try_olaps(reader,chr,win)
```
"""
function try_olaps(reader::BAM.Reader,chr::String,win::Array{Int64,1})::Array{BAM.Record,1}

    # NOTE: known bug that has to do w/ BioAlignments when close to end?
    records_olap =
    try
        collect(eachoverlap(reader, chr, win[1]:win[2]))
    catch
        Array{BAM.Record,1}()
    end

    return records_olap

end # try_olaps
"""
    `read_bam(BAM_PATH, CHR, FEAT_ST, FEAT_END, CPG_POS, PE)`

Function that reads in BAM file in `BAM_PATH` and returns methylation vectors
for those records that overlap with (1-based) genomic coordinates
chr:chrSt-chrEnd at cpg_pos. The information was taken from:

  XS: meth calls (https://github.com/FelixKrueger/Bismark/tree/master/Docs)
  XS: uniqueness (http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)

# Examples
```julia-repl
julia> read_bam(BAM_PATH,"chr1",30,80,[40,60],false)
```
"""
function read_bam(bam_path::String, chr::String, f_st::Int64, f_end::Int64,
                  cpg_pos::Array{Int64,1},chr_size::Int64,pe::Bool)::Array{Array{Int64,1},1}

    # Number of CpG sites is determined by that in the region
    N = length(cpg_pos)

    # Expand window in PE case to include pair even if it's outside the window
    exp_win = pe ? [max(1,f_st-75),min(chr_size,f_end+75)] : [f_st,f_end]

    # Get records overlapping window.
    reader = open(BAM.Reader, bam_path, index=bam_path*".bai")
    records_olap = try_olaps(reader,chr,exp_win)
    length(records_olap)>0 || return Array{Int64,1}[]

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
            meth_call = R1_call * "."^max(0,dist_st-length(R1_call))
            meth_call *= R2_call[max(1,(length(R1_call)+1-dist_st)):end]
            record.strand in ["OT","CTOT"] ? OFFSET=1 : OFFSET=2
            OFFSET -= BAM.position(record.R1)
        end

        # Cross positions of CpG sites if template contains CpG sites
        obs_cpgs = [x.offset for x in eachmatch(r"[zZ]",meth_call)].-OFFSET
        length(obs_cpgs)>0 || continue

        # If overlapping CpG sites then store and add to xobs
        olap_cpgs = findall(x-> x in obs_cpgs, cpg_pos)
        length(olap_cpgs)>0 || continue
        x = zeros(Int64,N)
        obs_cpgs = obs_cpgs[findall(x-> x in cpg_pos,obs_cpgs)] .+ OFFSET
        x[olap_cpgs] = reduce(replace,["Z"=>1,"z"=>-1],init=split(meth_call[obs_cpgs],""))
        push!(xobs,x)
            # elseif (maximum(obs_cpgs) >= minimum(cpg_pos)) ⊻
            #     (maximum(cpg_pos)>=minimum(obs_cpgs))
            #     # CpG site positions are not intercalated
            # else
            #     println(stderr,"[$(now())]: Read/s $(BAM.tempname(record.R1))" *
            #       " with flag $(BAM.flag(record.R1)) has unused CpG sites.")
            # end
            # end
    end # end loop over templates sequenced

    # Return
    return xobs

end # end read_bam
"""
    `write_gff!(OUT_GFF, GFF_RECORDS)`

Function that appends `GFF_RECORDS` into `OUT_GFF`.

# Examples
```julia-repl
julia> write_gff!(out_gff_path, gff_records)
```
"""
function write_gff!(out_gff_path::String, gff_records::Vector{GFF3.Record})

    # Append gff records
    output = open(out_gff_path, "a")
    writer = GFF3.Writer(output)
    while length(gff_records)>0
        write(writer,popfirst!(gff_records))
    end
    close(writer)

end # end write_gff
"""
    `next_record(READER, CHR_NAMES, RECORD)`

Recursive function that returns next record in stream if in chr_names, and false
if the stream reached its end. This function is only used in `read_vcf()`.

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
    `is_het_cpg!(VAR,SEQ,H1,H2)`

Function that adds position of heterozygous CpG site, if existant, to the
correct haplotype. Assumes that ref is the one pertaining to genome 1, and alt
is the one pertaining to genome 2.

# Examples
```julia-repl
julia> is_het_cpg!(var,seq,h1,h2)
```
"""
function is_het_cpg!(var::VCF.Record,seq::FASTA.Record,h1::Array{Int64,1},
                     h2::Array{Int64,1})
    # Initialize
    ref_var = join(VCF.ref(var))
    alt_var = join(VCF.alt(var))

    # Check if heterozygous CpG site in XG context
    if ((ref_var=="C") && (alt_var in ["A","G"]))
        # If followed by G, then heterozygous CpG in haplotype 1
        if FASTA.sequence(seq,VCF.pos(var):(VCF.pos(var)+1))[2]==DNA_G
            push!(h1,VCF.pos(var))
        end
    elseif ((alt_var=="C") && (ref_var in ["A","G"]))
        # If followed by G, then heterozygous CpG in haplotype 2
        if FASTA.sequence(seq,VCF.pos(var):(VCF.pos(var)+1))[2]==DNA_G
            push!(h2,VCF.pos(var))
        end
    end

    # Check if het CpG site in CX context
    if ((ref_var=="G") && (alt_var in ["A","C","T"]))
        # If preceeded by C, then heterozygous CpG in haplotype 1
        if FASTA.sequence(seq,(VCF.pos(var)-1):VCF.pos(var))[1]==DNA_C
            push!(h1,VCF.pos(var)-1)
        end
    elseif ((alt_var=="G") && (ref_var in ["A","C","T"]))
        # If preceeded by C, then heterozygous CpG in haplotype 2
        if FASTA.sequence(seq,(VCF.pos(var)-1):VCF.pos(var))[1]==DNA_C
            push!(h2,VCF.pos(var)-1)
        end
    end

end # end is_het_cpg!
"""
    `get_records_ps!(READER_VCF,SEQ,RECORD,CURR_PS,WIN_END,H1,H2)`

Recursive function that returns the end of current PS phased SNPs along the
next SNP. This function is only used in `read_vcf()`.

# Examples
```julia-repl
julia> get_records_ps!(reader_VCF,seq,record,curr_ps,wEnd,h1,h2)
```
"""
function get_records_ps!(reader_vcf::VCF.Reader,seq::FASTA.Record,
                         record::VCF.Record,curr_ps::String,wEnd::Int64,
                         h1::Array{Int64,1},h2::Array{Int64,1})::Tuple{VCF.Record,Int64}

    # Check if record has same PS
    ind_PS = findfirst(x->x=="PS",VCF.format(record))
    VCF.genotype(record)[1][ind_PS]==curr_ps || return record,wEnd

    # Check GT string and update het. CpG sites
    ind_GT = findfirst(x->x=="GT",VCF.format(record))
    curr_gt = VCF.genotype(record)[1][ind_GT]
    curr_gt in ["0/1","0|1"] ? is_het_cpg!(record,seq,h1,h2) : is_het_cpg!(record,seq,h2,h1)

    # Continue reading
    wEnd = VCF.pos(record)
    if !eof(reader_vcf)
        read!(reader_vcf, record)
        if "PS" in VCF.format(record)
            record,wEnd=get_records_ps!(reader_vcf,seq,record,curr_ps,wEnd,h1,h2)
        else
            return record,wEnd
        end
    else
        # return empty record if finished stream
        return VCF.Record(),wEnd
    end

end # end get_records_ps
"""
    `read_vcf(OUT_GFF_PATH,FASTA_PATH,VCF_PATH,WINDOW_SIZE)`

Function that creates a GFF file containing the heterozygous SNPs along with
the positions of the sorrounding CpG sites within a window of `WINDOW_SIZE`. For
phased VCF records, GT must be specified using A|B notation, while for unphased
records it assumes GT is always "0/1". The phasing is specified through the
standard notation PS.

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
    read!(reader_vcf, record)

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
    while true

        # Check if chromosome is in FASTA file. If not, obtain the next valid 1
        record = next_record(reader_vcf, chr_names, record)

        # Check we have loaded the right chromosome from FASTA
        if !(curr_chr == VCF.chrom(record))
            curr_chr = VCF.chrom(record)
            reader_fasta = open(FASTA.Reader,fasta_path,index=fasta_path*".fai")
            fasta_record = reader_fasta[curr_chr]
            chr_size = chr_sizes[findfirst(x->x==curr_chr, chr_names)]
            close(reader_fasta)
        end

        # Get entire (contiguous) haplotype using the PS field
        wSt = VCF.pos(record)
        wEnd = VCF.pos(record)
        het_cpgs1 = Array{Int64,1}()
        het_cpgs2 = Array{Int64,1}()
        if "PS" in VCF.format(record)
            # If PS tag, then phased
            ind_PS = findfirst(x->x=="PS",VCF.format(record))
            curr_ps = VCF.genotype(record)[1][ind_PS]
            record,wEnd = get_records_ps!(reader_vcf, fasta_record, record, curr_ps,
                                          wEnd, het_cpgs1, het_cpgs2)
        else
            # If no PS tag, then single SNP (marked as NOPS)
            curr_ps = "NOPS"
            is_het_cpg!(record,fasta_record,het_cpgs1,het_cpgs2)
            eof(reader_vcf) ? record=VCF.Record() : read!(reader_vcf, record)
        end

        # Obtain DNA sequence from reference and CpG loci from it
        wSt -= Int64(window_size/2)
        wEnd += Int64(window_size/2)
        win = [max(1, wSt), min(chr_size, wEnd)]
        wSeq = convert(String,FASTA.sequence(fasta_record,win[1]:win[2]))
        cpg_pos = map(x->getfield(x,:offset),eachmatch(r"CG",wSeq)).+win[1].-1

        # Store haplotype & corresponding CpG sites if homozygous CpG site/s
        if length(cpg_pos)>0
            out_str = "$(curr_chr)\t.\t$(curr_ps)\t$(win[1])\t$(win[2])"*
                      "\t.\t.\t.\tN=$(length(cpg_pos));CpGs=$(cpg_pos)"
            length(het_cpgs1)>0 && (out_str*=";hetCpGg1=$(het_cpgs1)")
            length(het_cpgs2)>0 && (out_str*=";hetCpGg2=$(het_cpgs2)")
            push!(gff_records,GFF3.Record(out_str))
        end

        # Check if gff_records object is too big (8 bytes x record) and dump it
        sizeof(gff_records)>=80000 && write_gff!(out_gff_path, gff_records)

        # If new record empty break while loop
        VCF.haspos(record) || break

    end

    # Dump remaining records
    write_gff!(out_gff_path, gff_records)

end # end read_vcf
"""
    `read_gff_chr(GFF_PATH,CHR)`

Function that reads in a GFF3 file in GFF_PATH and returns an array of
GenomicFeatures Record objects contained in chromosome CHR. The format is the
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
    `write_bG(BEDGRAPH_PATH, BEDGRAPH_RECORDS)`

Function that writes records in `BEDGRAPH_RECORDS` into `BEDGRAPH_PATH`.

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
    `get_cpg_pos(FEAT_ATTS)`

Function that takes the attributes in the GFF file and returns the homozygous
and heterozygous CpG sites in it.

# Examples
```julia-repl
julia> get_cpg_pos(f_atts)
```
"""
function get_cpg_pos(atts::Dict{String,Array{String,1}})::Array{Array{Int64,1},1}

    # Retrieve homozygous CpG sites
    hom = [parse(Int64,replace(s,r"[\[\] ]"=>"")) for s in atts[:"CpGs"]]

    # Retrieve heterozygous CpG sites for each haplotype (if existant)
    het1=Array{Int64,1}()
    het2=Array{Int64,1}()
    if haskey(atts,"hetCpGg1")
        het1=[parse(Int64,replace(s,r"[\[\] ]"=>"")) for s in atts[:"hetCpGg1"]]
    end
    if haskey(atts,"hetCpGg2")
        het2=[parse(Int64,replace(s,r"[\[\] ]"=>"")) for s in atts[:"hetCpGg2"]]
    end

    return [hom,sort!(append!(het1,hom)),sort!(append!(het2,hom))]

end # end get_cpg_pos
"""
    `get_ns(CPG_POS,BLOCK_SIZE)`

Functions that returns [N1,...,NK] given the position of the CpG sites and
a block size BLOCK_SIZE.

# Examples
```julia-repl
julia> get_ns([100,250,300,350],200)
2-element Array{Int64,1}:
 3
 1
```
"""
function get_ns(cpg_pos::Array{Int64,1},b_size::Int64)::Array{Int64,1}

    # Partition array
    y = Int64[]
    for p in cpg_pos
        j = 0
        while true
            ((p-cpg_pos[1]-j*b_size)>=0) && ((p-cpg_pos[1]-(j+1)*b_size)<=0) && break
            j += 1
        end
        push!(y,j)
    end

    # Return array of arrays
    return [sum(y.==i) for i in unique(y)]

end # end get_ns
"""
    `mean_cov(XOBS)`

Function returns the average coverage per CpG given some observations.

# Examples
```julia-repl
julia> mean_cov(xobs)
```
"""
function mean_cov(xobs::Array{Array{Int64,1},1})::Float64

    # Return 0 if no observations
    length(xobs)>0 || return 0.0
    return sum(sum(map(x->abs.(x),xobs),dims=(1))[1])/length(xobs[1])

end # end mean_cov
"""
    `run_asm_analysis(BAM1_PATH,BAM2_PATH,VCF_PATH,FASTA_PATH,WINDOW_SIZE,OUT_PATH)`

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
                          pe::Bool=true,b_size::Int64=100)

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
    mml1_path = "$(out_path)/$(prefix_sample)_mml1.bedGraph"
    mml2_path = "$(out_path)/$(prefix_sample)_mml2.bedGraph"
    h1_path = "$(out_path)/$(prefix_sample)_h1.bedGraph"
    h2_path = "$(out_path)/$(prefix_sample)_h2.bedGraph"
    mi_path = "$(out_path)/$(prefix_sample)_mi.bedGraph"
    pVal_path = "$(out_path)/$(prefix_sample)_pVal.bedGraph"

    # Check for existance of at least an output files
    if isfile.(mml1_path,mml2_path,h1_path,h2_path,mi_path,pVal_path)
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
    h1_recs = Array{Tuple{String,Int64,Int64,Float64},1}()
    h2_recs = Array{Tuple{String,Int64,Int64,Float64},1}()
    mi_recs = Array{Tuple{String,Int64,Int64,Float64},1}()
    mml1_recs = Array{Tuple{String,Int64,Int64,Float64},1}()
    mml2_recs = Array{Tuple{String,Int64,Int64,Float64},1}()
    pVal_recs = Array{Tuple{String,Int64,Int64,Float64},1}()

    # Loop over chromosomes
    n = n1 = n2 = []
    ex = exx = mi = 0.0
    cpg_pos = Array{Array{Int64,1},1}()
    f_atts = Dict{String,Array{String,1}}()
    tot_feats = int_feats_1 = int_feats_2 = int_feats_mi = 0
    for chr in chr_names

        # Get windows pertaining to current chromosome
        println(stderr,"[$(now())]: Processing 🧬  $(chr) ...")
        features_chr = read_gff_chr(out_gff_path,chr)
        chr_size = chr_sizes[findfirst(x->x==chr, chr_names)]

        # Loop over windows in chromosome
        @showprogress 1 "Computing..." for feat in features_chr
            # Get window of interest
            tot_feats += 1
            f_st = GFF3.seqstart(feat)
            f_end = GFF3.seqend(feat)
            f_atts = Dict(GFF3.attributes(feat))

            # Get homozygous & heterozygous CpG sites (1:hom, 2:hap1, 3: hap2)
            cpg_pos = get_cpg_pos(f_atts)

            # Get vectors from BAM1/2 overlapping feature
            n = get_ns(cpg_pos[1],b_size)
            n1 = get_ns(cpg_pos[2],b_size)
            n2 = get_ns(cpg_pos[3],b_size)
            xobs1 = read_bam(bam1_path,chr,f_st,f_end,cpg_pos[2],chr_size,pe)
            xobs2 = read_bam(bam2_path,chr,f_st,f_end,cpg_pos[3],chr_size,pe)

            # Estimate each single-allele model, mml and h
            if mean_cov(xobs1)>=THRESH_COV
                eta1 = est_eta(n1,xobs1)
                ex = comp_ex(n1,eta1[1:length(n1)],eta1[end])
                exx = comp_exx(n1,eta1[1:length(n1)],eta1[end])
                h1 = comp_shanH(n1,eta1[1:length(n1)],eta1[end],ex,exx)
                push!(mml1_recs,(chr,f_st,f_end,comp_mml(ex)))
                push!(h1_recs,(chr,f_st,f_end,h1))
                int_feats_1 += 1
            end
            if mean_cov(xobs2)>=THRESH_COV
                eta2 = est_eta(n2,xobs2)
                ex = comp_ex(n2,eta2[1:length(n2)],eta2[end])
                exx = comp_exx(n2,eta2[1:length(n2)],eta2[end])
                h2 = comp_shanH(n2,eta2[1:length(n2)],eta2[end],ex,exx)
                push!(mml2_recs,(chr,f_st,f_end,comp_mml(ex)))
                push!(h2_recs,(chr,f_st,f_end,h2))
                int_feats_2 += 1
            end

            # Compute mutual information
            if (mean_cov(xobs1)>=THRESH_COV) && (mean_cov(xobs2)>=THRESH_COV)
                ind1 = findall(x->x==true,[p in cpg_pos[1] for p in cpg_pos[2]])
                ind2 = findall(x->x==true,[p in cpg_pos[1] for p in cpg_pos[3]])
                xobs1 = [x[ind1] for x in xobs1]
                xobs2 = [x[ind2] for x in xobs2]
                eta = est_eta(n,vcat(xobs1,xobs2))
                ex = comp_ex(n,eta[1:length(n)],eta[end])
                exx = comp_exx(n,eta[1:length(n)],eta[end])
                mi = comp_mi(comp_shanH(n,eta[1:length(n)],eta[end],ex,exx),h1,h2)
                push!(mi_recs,(chr,f_st,f_end,mi))
                # push!(pVal_recs,(chr,f_st,f_end,perm_test(xobs1,xobs2,mi,cpg_pos)))
                int_feats_mi += 1
            end

        end

        # Add values to respective output BigWig files after each chr
        mml1_recs = write_bG(mml1_recs, mml1_path)
        mml2_recs = write_bG(mml2_recs, mml2_path)
        h1_recs = write_bG(h1_recs, h1_path)
        h2_recs = write_bG(h2_recs, h2_path)
        mi_recs = write_bG(mi_recs, mi_path)
        pVal_recs = write_bG(pVal_recs, pVal_path)

    end

    # Print number of features interrogated and finished message
    println(stderr,"[$(now())]: Total number of features: $(tot_feats).")
    println(stderr,"[$(now())]: G1:$(round(100*int_feats_1/tot_feats; sigdigits=2))%.")
    println(stderr,"[$(now())]: G2:$(round(100*int_feats_2/tot_feats;sigdigits=2))%.")
    println(stderr,"[$(now())]: Both:$(round(100*int_feats_mi/tot_feats;sigdigits=2))%.")
    print(stderr,"[$(now())]: Done. 😄")

end # end run_asm_analysis
