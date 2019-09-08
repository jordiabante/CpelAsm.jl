###################################################################################################
# CONSTANTS
###################################################################################################
const BG_BUFFER = 50000                     # Bytes of bedGraph records until write
const GFF_BUFFER = 500000                   # Bytes of GFF records until write
const THRESH_MAPQ = 39                      # MAPQ threshold (-10*log(p)) only true uni-reads
const FLAGS_ALLOWED = [0,16,83,99,147,163]  # Flags allowed in BAM recs
###################################################################################################
# STRUCTS
###################################################################################################
struct AlignTemp
    strand::String      # Methylation call strand
    R1::BAM.Record      # First record from left to right
    R2::BAM.Record      # Second record from left to right
end
mutable struct AllAlignTemps
    paired::Bool
    templates::Array{AlignTemp,1}
end
###################################################################################################
# FUNCTIONS
###################################################################################################
"""
    `get_align_strand(PAIRED_END,FLAG1,FLAG2)`

Function that returns the strand of methylation call based on Bismark's logic. In single-end mode,
OT and CTOT reads will both receive a FLAG of 0 while OB and CTOB both receive a FLAG of 16. In
paired end mode:

             Read 1       Read 2
  OT:         99          147
  OB:         83          163
  CTOT:      147           99
  CTOB:      163           83

This table was extracted from

    https://github.com/FelixKrueger/Bismark/issues/151

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
            println(stderr,"[$(now())]: Exiting julia ...")
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
            println(stderr, "[$(now())]: Exiting julia ...")
            exit(1)
        end
    end

    # Return strand
    return s

end # end get_align_strand
"""
    `order_bams(PAIRED_END,RECORDS)`

Function that returns an AlignTemp object with R1 as the first record in the forward strand and
with the methylation call strand taken from `get_align_strand`.

# Examples
```julia-repl
julia> JuliASM.order_bams(true,RECORDS)
```
"""
function order_bams(pe::Bool,records::Vector{BAM.Record})::AlignTemp

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

Function that takes a set of records and returns an AllAlignTemps object that contains all the
properly aligned templates as an array of AlignTemp, which contains information about the
methylation call strand as well as the relevant BAM records. In the PE case, R1 corresponds to
the BAM record that appears before in the forward strand.

# Examples
```julia-repl
julia> JuliASM.clean_records(true,RECORDS)
```
"""
function clean_records(pe::Bool,records::Vector{BAM.Record})::AllAlignTemps

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
        end
    else
        # Handle SE case
        out.templates = [AlignTemp(get_align_strand(pe,BAM.flag(x),UInt16(0)),x,BAM.Record()) for
            x in records]
    end

    # Return struct
    return out

end # end clean_records
"""
    `try_olaps(READER,CHR,WINDOW)`

Function that tries to find BAM records overlaping with `CHR` at positions `WINDOW`.

# Examples
```julia-repl
julia> JuliASM.try_olaps(reader,chr,win)
```
"""
function try_olaps(reader::BAM.Reader,chr::String,win::Vector{Int64})::Vector{BAM.Record}

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
    `read_bam(BAM_PATH,CHR,FEAT_ST,FEAT_END,CPG_POS,PE,TRIM)`

Function that reads in BAM file in `BAM_PATH` and returns methylation vectors for those records
that overlap with (1-based) genomic coordinates `chr:chrSt-chrEnd` at `cpg_pos`. A trimming given
by TRIM=(Rf_5p,Rf_3p,Rr_5p,Rr_3p) is applied to the reads. The information was taken from:

  XS: meth calls (https://github.com/FelixKrueger/Bismark/tree/master/Docs)
  XS: uniqueness (http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)

For info on OT, CTOT, OB, CTOB nomenclature see

    http://software.broadinstitute.org/software/igv/book/export/html/37.

# Examples
```julia-repl
julia> JuliASM.read_bam(BAM_PATH,"chr1",30,80,[40,60],false,(0,0,0,0))
```
"""
function read_bam(bam::String,chr::String,hap_st::Int64,hap_end::Int64,cpg_pos::Vector{Int64},
                  chr_size::Int64,pe::Bool,trim::NTuple{4,Int64})::Array{Vector{Int64},1}

    # Number of CpG sites is determined by that in the region
    N = length(cpg_pos)

    # Expand window in PE case to include pair even if it's outside the window
    exp_win = pe ? [max(1,hap_st-75),min(chr_size,hap_end+75)] : [hap_st,hap_end]

    # Get records overlapping window.
    reader = open(BAM.Reader,bam,index=bam*".bai")
    records_olap = try_olaps(reader,chr,exp_win)
    length(records_olap)>0 || return Vector{Int64}[]
    close(reader)

    # Relevant flags in BAM file (both Bismark and Arioc)
    filter!(x-> (BAM.ismapped(x)) && (haskey(x,"XM")) && (!haskey(x,"XS")) && (BAM.flag(x) in
           FLAGS_ALLOWED) && (BAM.mappingquality(x)>THRESH_MAPQ),records_olap)

    # Clean records & organize them
    records_org = clean_records(pe,records_olap)

    # Loop over organized records
    xobs = Vector{Int64}[]
    for record in records_org.templates
        # Get methylation call and offset (depends on strand where call is made)
        if !(records_org.paired)
            # Obtain methylation call from single end ()
            meth_call = record.R1[:"XM"]
            meth_call = SubString(meth_call,(1+trim[1]):(length(meth_call)-trim[2]))
            OFFSET= record.strand in ["OT","CTOT"] ? 1 : 2
            OFFSET -= BAM.position(record.R1)+trim[1]
        else
            # Obtain methylation call
            R1_call = record.R1[:"XM"]
            R2_call = record.R2[:"XM"]
            R1_call = SubString(R1_call,(1+trim[1]):(length(R1_call)-trim[2]))
            R2_call = SubString(R2_call,(1+trim[4]):(length(R2_call)-trim[3]))
            dist_st = abs(BAM.position(record.R2)+trim[4] - BAM.position(record.R1)-trim[1])
            meth_call = R1_call * "."^max(0,dist_st-length(R1_call))
            meth_call *= R2_call[max(1,(length(R1_call)+1-dist_st)):end]
            OFFSET= record.strand in ["OT","CTOT"] ? 1 : 2
            OFFSET -= BAM.position(record.R1)+trim[1]
        end

        # Cross positions of CpG sites if template contains CpG sites
        obs_cpgs = [x.offset for x in eachmatch(r"[zZ]",meth_call)] .- OFFSET
        length(obs_cpgs)>0 || continue

        # If overlapping CpG sites then store and add to xobs
        olap_cpgs = findall(x-> x in obs_cpgs, cpg_pos)
        length(olap_cpgs)>0 || continue
        x = zeros(Int64,N)
        obs_cpgs = obs_cpgs[findall(x-> x in cpg_pos,obs_cpgs)] .+ OFFSET
        x[olap_cpgs] = reduce(replace,["Z"=>1,"z"=>-1],init=split(meth_call[obs_cpgs],""))
        push!(xobs,x)

    end # end loop over templates sequenced

    # Return
    return xobs

end # end read_bam
"""
    `write_gff!(OUT_GFF,GFF_RECORDS)`

Function that appends `GFF_RECORDS` into `OUT_GFF`.

# Examples
```julia-repl
julia> JuliASM.write_gff!(gff, gff_records)
```
"""
function write_gff!(gff::String,gff_records::Vector{GFF3.Record})

    # Append gff records
    output = open(gff, "a")
    writer = GFF3.Writer(output)
    while length(gff_records)>0
        write(writer,popfirst!(gff_records))
    end
    close(writer)

    # Return
    return nothing

end # end write_gff
"""
    `next_record(READER,CHR_NAMES,RECORD)`

Recursive function that returns next record in stream if in chr_names, and false if the stream
reached its end. This function is only used in `gen_gffs()`.

# Examples
```julia-repl
julia> JuliASM.next_record(reader, chr_names, record)
```
"""
function next_record(reader::VCF.Reader,chr_names::Vector{String},record::VCF.Record)::VCF.Record

    # Check if record in chr_names and end of reader
    if VCF.chrom(record) in chr_names
        # If we're good, return record
        return record
    elseif !eof(reader)
        # If not end of file, call function again
        read!(reader, record)
        next_record(reader, chr_names, record)
    else
        # return empty record if finished stream
        return VCF.Record()
    end

end # end next_record
"""
    `is_het_cpg!(VAR,SEQ,H1,H2)`

Function that adds position of heterozygous CpG site, if existant, to the correct haplotype.
Assumes that REF is the one pertaining to genome 1, and ALT is the one pertaining to genome 2.

# Examples
```julia-repl
julia> JuliASM.is_het_cpg!(var,seq,h1,h2)
```
"""
function is_het_cpg!(var::VCF.Record,seq::FASTA.Record,h1::Vector{Int64},h2::Vector{Int64})
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

    # Return
    return nothing

end # end is_het_cpg!
"""
    `get_records_ps!(VCF_READER,SEQ,RECORD,CURR_PS,WIN_END,H1,H2)`

Recursive function that returns the end of current PS phased SNPs along the next SNP. This
function is only used in `gen_gffs()`.

# Examples
```julia-repl
julia> JuliASM.get_records_ps!(reader,seq,record,curr_ps,wEnd,h1,h2)
```
"""
function get_records_ps!(reader::VCF.Reader,seq::FASTA.Record,record::VCF.Record,curr_ps::String,
                         wEnd::Int64,h1::Vector{Int64},h2::Vector{Int64})::Tuple{VCF.Record,Int64}

    # Check if record has same PS
    ind_PS = findfirst(x->x=="PS",VCF.format(record))
    VCF.genotype(record)[1][ind_PS]==curr_ps || return record,wEnd

    # Check GT string and update het. CpG sites
    ind_GT = findfirst(x->x=="GT",VCF.format(record))
    curr_gt = VCF.genotype(record)[1][ind_GT]
    curr_gt in ["0/1","0|1"] ? is_het_cpg!(record,seq,h1,h2) : is_het_cpg!(record,seq,h2,h1)

    # Continue reading
    wEnd = VCF.pos(record)
    if !eof(reader)
        read!(reader, record)
        if "PS" in VCF.format(record)
            record,wEnd = get_records_ps!(reader,seq,record,curr_ps,wEnd,h1,h2)
        else
            return record,wEnd
        end
    else
        # return empty record if finished stream
        return VCF.Record(),wEnd
    end

end # end get_records_ps
"""
    `expand_win(FIRST_SNP,LAST_SNP,WIN_EXP,BLOCK_SIZE,CHR_SIZE)`

Function that returns an expanded window given the position of the first and last SNPs, the
expansion to be done to the left and right, of the first and last SNP respectively, to include
nearby CpG sites, and the size of the subregions into which the whole region will be divided into.

# Examples
```julia-repl
julia> JuliASM.expand_win(50,80,10,2000)
2-element Array{Int64,1}:
 40
 90
```
"""
function expand_win(l_snp::Int64,r_snp::Int64,win_exp::Int64,chr_size::Int64)::Vector{Int64}

    # Return extended window
    return [max(1,l_snp-win_exp),min(chr_size,r_snp+win_exp)]

end # end expand_win
"""
    `find_num_models(CPG_POS,NMAX,NMOD)`

Function that returns the number of models required so as to avoid joint model with N>NMAX.

# Examples
```julia-repl
julia> JuliASM.find_num_models([50,55,60,65],2,1)
2
```
"""
function find_num_models(cpg_pos::Vector{Int64},n_max::Int64,n_mod::Int64)::Int64

    # Check if resulting models are acceptable
    return length(cpg_pos)/n_mod>n_max ? find_num_models(cpg_pos,n_max,n_mod+1) : n_mod

end # end find_num_models
"""
    `gen_gffs([HET_GFF_PATH,HOM_GFF_PATH],FA_PATH,VCF_PATH,WIN_EXP,NMAX)`

Function that creates a GFF file containing the genomic regions to be analyzed by `JuliASM`. If a
set of SNPs is phased, this function will create a single window for all. In case, the set of SNPs
is not phased, then an individual window will be created for each SNP. The phasing of the VCF
records should be specified by means of the standard `PS` tag. In addition, a GFF file containing
the homozygous portion of the DNA is also generated together with the position of the CpG sites
in it.

# Examples
```julia-repl
julia> JuliASM.gen_gffs([HET_GFF_PATH,HOM_GFF_PATH],FA_PATH,VCF_PATH,WIN_EXP,NMAX)
```
"""
function gen_gffs(gff::Vector{String},fa::String,vcf::String,win_exp::Int64,n_max::Int64)

    # Initialize VCF variables
    het_records = Vector{GFF3.Record}()
    hom_records = Vector{GFF3.Record}()
    reader_vcf = open(VCF.Reader,vcf)
    record = VCF.Record()
    read!(reader_vcf,record)

    # Initialize FASTA variables
    prev_end = 0
    reader_fa = open(FASTA.Reader, fa, index=fa*".fai")
    chr_names = reader_fa.index.names
    curr_chr = chr_names[1]
    fa_record = reader_fa[curr_chr]
    close(reader_fa)

    # Chrom sizes
    chr_sizes = reader_fa.index.lengths
    chr_size = chr_sizes[findfirst(x->x==curr_chr,chr_names)]

    # Loop over variants
    while true

        # Check if chromosome is in FASTA file. If not, obtain the next valid 1
        record = next_record(reader_vcf,chr_names,record)

        # Check we have loaded the right chromosome from FASTA
        if !(curr_chr == VCF.chrom(record))
            curr_chr = VCF.chrom(record)
            reader_fa = open(FASTA.Reader,fa,index=fa*".fai")
            fa_record = reader_fa[curr_chr]
            chr_size = chr_sizes[findfirst(x->x==curr_chr,chr_names)]
            close(reader_fa)
            prev_end=0
        end

        # Get entire (contiguous) haplotype using the PS field
        wSt = VCF.pos(record)
        wEnd = VCF.pos(record)
        het1 = Vector{Int64}()
        het2 = Vector{Int64}()
        if "PS" in VCF.format(record)
            # If PS tag, then phased
            ind_PS = findfirst(x->x=="PS",VCF.format(record))
            curr_ps = VCF.genotype(record)[1][ind_PS]
            record,wEnd = get_records_ps!(reader_vcf,fa_record,record,curr_ps,wEnd,het1,het2)
        else
            # If no PS tag, then single SNP (marked as NOPS)
            curr_ps = "NOPS"
            is_het_cpg!(record,fa_record,het1,het2)
            eof(reader_vcf) ? record=VCF.Record() : read!(reader_vcf, record)
        end

        # Expand window win_exp on each side
        win = expand_win(wSt,wEnd,win_exp,chr_size)

        # Obtain DNA sequence from reference and homozygous CpG loci from it
        wSeq = convert(String,FASTA.sequence(fa_record,win[1]:win[2]))
        cpg_pos = map(x->getfield(x,:offset),eachmatch(r"CG",wSeq)).+win[1].-1

        # Store heterozygous haplotype & corresponding CpG sites if homozygous CpG site/s
        if length(cpg_pos)>0

            # Obtain required number of models based on homozygous CpG sites
            n_mod = find_num_models(cpg_pos,n_max,1)
            mod_size = div(length(cpg_pos),n_mod)
            mod_rem = length(cpg_pos)-n_mod*mod_size

            # Print each set separately
            for i=1:n_mod
                # Index of subset of CpG sites in cpg_pos
                hom_rng = (mod_size*(i-1)+1):min((mod_size*i),length(cpg_pos))
                # Add remainders to last block
                i==n_mod && mod_rem>0 && (hom_rng=extrema(hom_rng)[1]:length(cpg_pos))
                # Get subset of CpG sites for each
                hom_mod = cpg_pos[hom_rng]
                # Get submodel boundaries
                inf_bound = i==1 ? win[1] : cpg_pos[max(1,extrema(hom_rng)[1]-1)]
                sup_bound = i==n_mod ? win[2] : hom_mod[end]
                # Get heterozygous CpG sites in current submodel
                het1_mod = het1[findall(x-> inf_bound <= x <= sup_bound,het1)]
                het2_mod = het2[findall(x-> inf_bound <= x <= sup_bound,het2)]
                # Make submodel window a multiple of G (no new CpGs included)
                win_mod = expand_win(inf_bound,sup_bound,0,chr_size)
                # Push output string
                out_str = "$(curr_chr)\t.\t$(curr_ps)-$(i)/$(n_mod)\t$(win_mod[1])\t$(win_mod[2])"
                out_str *= "\t.\t.\t.\tN=$(length(hom_mod));CpGs=$(hom_mod)"
                length(het1_mod)>0 && (out_str*=";hetCpGg1=$(het1_mod)")
                length(het2_mod)>0 && (out_str*=";hetCpGg2=$(het2_mod)")
                # Push GFF record
                push!(het_records,GFF3.Record(out_str))
            end

        end

        # Check if het_records object is too big and empty it if so
        sizeof(het_records)>GFF_BUFFER && write_gff!(gff[1],het_records)

        # Obtain homozygous portion
        hom_win = [prev_end+1,win[1]-1]
        hom_win = [max(1,hom_win[1]),min(chr_size,hom_win[2])]

        # Obtain DNA sequence from reference and homozygous CpG loci from it
        wSeq = convert(String,FASTA.sequence(fa_record,hom_win[1]:hom_win[2]))
        cpg_pos = map(x->getfield(x,:offset),eachmatch(r"CG",wSeq)).+hom_win[1].-1

        # Store homozygous stretch of CpG sites found
        if length(cpg_pos)>0
            out_str = "$(curr_chr)\t.\t.\t$(hom_win[1])\t$(hom_win[2])\t$(length(cpg_pos))\t.\t.\t"
            out_str *= "CpGs=$(cpg_pos)"
            push!(hom_records,GFF3.Record(out_str))
        end

        # Update previous end
        prev_end = win[2]

        # Check if hom_records object is too big and empty it if so
        sizeof(hom_records)>GFF_BUFFER && write_gff!(gff[2],hom_records)

        # If new record empty break while loop
        VCF.haspos(record) || break

    end

    # Dump remaining records
    write_gff!(gff[1],het_records)
    write_gff!(gff[2],hom_records)

    # Return
    return nothing

end # end gen_gffs
"""
    `read_gff_chr(GFF_PATH,CHR)`

Function that reads in a GFF3 file in `GFF_PATH` and returns an array of GenomicFeatures Record
objects contained in chromosome `CHR`. The format is the standard defined by ENSEMBL in
https://useast.ensembl.org/info/website/upload/gff3.html.

# Examples
```julia-repl
julia> JuliASM.read_gff_chr(GFF_PATH,"chr1")
```
"""
function read_gff_chr(gff::String,chr::String)::Vector{GFF3.Record}

    # Load genomic features from a GFF3 file
    features = open(collect,GFF3.Reader,gff)

    # Keep features in chr
    filter!(x -> GFF3.seqid(x)==chr, features)

    # Return
    return features

end # end read_gff_chr
"""
    `write_tobs(RECORDS,CHR,PATH)`

Function that writes records in `RECORDS` into `PATH`.

# Examples
```julia-repl
julia> JuliASM.write_tobs(RECORDS,CHR,PATH)
```
"""
function write_tobs(recs::Vector{Tuple{Int64,Int64,Float64,Int64,Int64}},chr::String,path::String)

    # Open output bedGraph file in append mode (no need to close it)
    open(path, "a") do f
        for i=1:length(recs)
            recs[i][1]!=0 || continue
            write(f,"$(chr)\t"*join(string.(collect(recs[i])),"\t"),"\n")
        end
    end

    # Return
    return nothing

end # end write_tobs
"""
    `get_cpg_pos(FEAT_ATTS)`

Function that takes the attributes in the GFF file and returns the homozygous and heterozygous CpG
sites in it.

# Examples
```julia-repl
julia> JuliASM.get_cpg_pos(feat_atts)
```
"""
function get_cpg_pos(atts::Dict{String,Vector{String}})::Array{Vector{Int64},1}

    # Retrieve homozygous CpG sites
    hom = [parse(Int64,replace(s,r"[\[\] ]"=>"")) for s in atts[:"CpGs"]]

    # Retrieve heterozygous CpG sites for each haplotype (if existant)
    het1=Vector{Int64}()
    het2=Vector{Int64}()
    if haskey(atts,"hetCpGg1")
        het1=[parse(Int64,replace(s,r"[\[\] ]"=>"")) for s in atts[:"hetCpGg1"]]
    end
    if haskey(atts,"hetCpGg2")
        het2=[parse(Int64,replace(s,r"[\[\] ]"=>"")) for s in atts[:"hetCpGg2"]]
    end

    return [hom,sort!(append!(het1,hom)),sort!(append!(het2,hom))]

end # end get_cpg_pos
"""
    `get_ns(CPG_POS,G_MAX,HAP_ST,HAP_END)`

Functions that returns [N1,...,NK] given the position of the CpG sites & maximum granularity G_MAX.

# Examples
```julia-repl
julia> JuliASM.get_ns([100,200,300,350],200,90,400)
2-element Array{Int64,1}:
 2
 2
```
"""
function get_ns(cpg_pos::Vector{Int64},g_max::Int64,hap_st::Int64,hap_end::Int64)::Vector{Int64}

    # Check if need to partition
    (hap_end-hap_st)<=g_max && return [length(cpg_pos)]

    # Get K of model
    k = ceil((hap_end-hap_st+1)/g_max)

    # Get delimiters
    delimiters = vcat(collect(hap_st:(hap_end-hap_st+1)/k:hap_end),hap_end)

    # Partition CpG sites
    ids = [findfirst(y->y>=x,delimiters)-1 for x in cpg_pos]

    # Return [N1,...,NK] vector
    return [sum(ids.==i) for i in unique(ids)]

end # end get_ns
"""
    `mean_cov(XOBS)`

Function returns the average coverage per CpG given some observations `XOBS`.

# Examples
```julia-repl
julia> xobs=[[1,-1] for i=1:10]; append!(xobs,[[1,0] for i=1:10]);
julia> JuliASM.mean_cov(xobs)
15.0
```
"""
function mean_cov(xobs::Array{Vector{Int64},1})::Float64

    # Return 0 if no observations
    return length(xobs)>0 ? norm(hcat(xobs...),1)/length(xobs[1]) : 0.0

end # end mean_cov
"""
    `mean_cov_sr(XOBS)`

Function returns the average coverage per CpG given some observations `XOBS` per subregion.

# Examples
```julia-repl
julia> n=[1,1]; xobs=[[1,-1] for i=1:10]; append!(xobs,[[1,0] for i=1:10]);
julia> JuliASM.mean_cov_sr(xobs,n)
15.0
```
"""
function mean_cov_sr(xobs::Array{Vector{Int64},1},n::Vector{Int64})::Vector{Float64}

    # Obtain columns
    rows = [(sum(n[1:(i-1)])+1):sum(n[1:i]) for i=2:length(n)]
    pushfirst!(rows,1:n[1])

    # Return 0 if no observations
    return length(xobs)>0 ? [norm(hcat(xobs...)[rows[i],:],1)/n[i] for i=1:length(rows)] :
        zeros(Float64,length(n))

end # end mean_cov_sr
"""
    `get_outpaths(OUTDIR,PREFIX)`

Function that given outdir and prefix returns paths.

# Examples
```julia-repl
julia> JuliASM.get_outpaths(outdir,prefix)
```
"""
function get_outpaths(outdir::String,prefix::String)

    # Paths
    tobs_path = "$(outdir)/$(prefix)" .* ["_mml1","_mml2","_nme1","_nme2","_uc"] .* ".bedGraph"
    tnull_path = "$(outdir)/$(prefix)" .* ["_dmml_null","_dnme_null","_uc_null"] .* ".bedGraph"
    p_path = "$(outdir)/$(prefix)" .* ["_dmml_pvals","_dnme_pvals","_uc_pvals"] .* ".bedGraph"

    # Return path for Tobs and Tnull
    return tobs_path,tnull_path,p_path

end # end get_outpaths
"""
    `sort_bedgraphs(BEDGRAPH_FILES)`

Function that sorts bedGraph files in vector BEDGRAPH_FILES.

# Examples
```julia-repl
julia> JuliASM.sort_bedgraphs(bg_files)
```
"""
function sort_bedgraphs(bg_files::Vector{String})

    # Loop over files
    for f in bg_files
        run(`sort -V -k 1,1 -k 2,2n -o $f $f`)
    end

    # Return
    return nothing

end # end sort_bedgraphs
"""
    `print_log(MESSAGE)`

Function that prints MESSAGE to stderr.

# Examples
```julia-repl
julia> JuliASM.print_log("Hello")
Hello
```
"""
function print_log(mess::String)

    # Right format for date
    date = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println(stderr,"[$(date)]: "  * mess)
    flush(stderr)

    # Return
    return nothing

end # end print_log
"""
    `proc_obs_hap()`

Function that computes an observed dMML, dNME, and UC given a homozygous haplotype.

# Examples
```julia-repl
julia> proc_obs_hap()
```
"""
function proc_obs_hap(hap::GFF3.Record,chr::String,chr_size::Int64,bam1::String,bam2::String,
                      gff::String,fa::String,out_paths::Vector{String},pe::Bool,g_max::Int64,
                      cov_ths::Int64,trim::NTuple{4,Int64})::Vector{Tuple{Int64,Int64,Float64,
                      Int64,Int64}}

    # Empty output
    nan_out = [(0,0,0.0,0,0),(0,0,0.0,0,0),(0,0,0.0,0,0),(0,0,0.0,0,0),(0,0,0.0,0,0)]

    # Get window
    hap_st = GFF3.seqstart(hap)
    hap_end = GFF3.seqend(hap)

    # Get homozygous & heterozygous CpG sites (1:hom, 2:hap1, 3: hap2)
    cpg_pos = get_cpg_pos(Dict(GFF3.attributes(hap)))
    n = get_ns(cpg_pos[1],g_max,hap_st,hap_end)
    n1 = get_ns(cpg_pos[2],g_max,hap_st,hap_end)
    n2 = get_ns(cpg_pos[3],g_max,hap_st,hap_end)
    length(n)>0 || return nan_out

    # Get vectors from BAM1/2 overlapping haplotype & compute average coverage
    xobs1 = read_bam(bam1,chr,hap_st,hap_end,cpg_pos[2],chr_size,pe,trim)
    cov_ths <= mean_cov(xobs1) <= 400 || return nan_out
    xobs2 = read_bam(bam2,chr,hap_st,hap_end,cpg_pos[3],chr_size,pe,trim)
    cov_ths <= mean_cov(xobs2) <= 400 || return nan_out

    # Estimate each single-allele model and check if on boundary of parameter space
    θ1 = est_theta_sa(n1,xobs1)
    check_boundary(θ1) && return nan_out
    θ2 = est_theta_sa(n2,xobs2)
    check_boundary(θ2) && return nan_out

    # Get binary vector with homozygous CpG sites
    z1 = BitArray([p in cpg_pos[1] ? true : false for p in cpg_pos[2]])
    z2 = BitArray([p in cpg_pos[1] ? true : false for p in cpg_pos[3]])

    # Estimate intermediate quantities
    if all(z1)
        # In case all are homozygous CpG sites
        ∇1 = get_grad_logZ(n1,θ1)
    else
        # In case NOT all are homozygous CpG sites
        ex1 = comp_ex(n1,θ1[1:(end-1)],θ1[end])
        exx1 = comp_exx(n1,θ1[1:(end-1)],θ1[end])
    end
    if all(z2)
        # In case all are homozygous CpG sites
        ∇2 = get_grad_logZ(n2,θ2)
    else
        # In case NOT all are homozygous CpG sites
        ex2 = comp_ex(n2,θ2[1:(end-1)],θ2[end])
        exx2 = comp_exx(n2,θ2[1:(end-1)],θ2[end])
    end

    # Compute output
    mml1 = all(z1) ? comp_mml_∇(n1,∇1) : comp_mml(z1,ex1)
    mml2 = all(z2) ? comp_mml_∇(n2,∇2) : comp_mml(z2,ex2)
    nme1 = all(z1) ? comp_nme_∇(n1,θ1,∇1) : comp_nme(z1,n1,θ1[1:(end-1)],θ1[end],ex1,exx1)
    nme2 = all(z2) ? comp_nme_∇(n2,θ2,∇2) : comp_nme(z2,n2,θ2[1:(end-1)],θ2[end],ex2,exx2)
    uc = comp_uc(z1,z2,n1,n2,θ1,θ2,nme1,nme2)

    # Report positions based on CpG sites not hap_st, hap_end
    bed_st = minimum(minimum.(cpg_pos))
    bed_end = maximum(maximum.(cpg_pos))

    # Return output
    return [(bed_st,bed_end,mml1,sum(n),length(n1)),
            (bed_st,bed_end,mml2,sum(n),length(n2)),
            (bed_st,bed_end,nme1,sum(n),length(n1)),
            (bed_st,bed_end,nme2,sum(n),length(n2)),
            (bed_st,bed_end,uc,sum(n),length(n))]

end
"""
    `comp_tobs(BAM1_PATH,BAM2_PATH,GFF_PATH,FA_PATH,OUT_PATHS)`

Function that computes MML1/2, NME1/2, and UC on a pair of BAM files (allele 1, allele 2) given a
(phased) VCF file that contains the heterozygous SNPs and a FASTA file that contains the reference
genome.

# Examples
```julia-repl
julia> comp_tobs(BAM1_PATH,BAM2_PATH,GFF_PATH,FA_PATH,OUT_PATHS)
```
"""
function comp_tobs(bam1::String,bam2::String,gff::String,fa::String,out_paths::Vector{String};
                   pe::Bool=true,g_max::Int64=300,cov_ths::Int64=5,trim::NTuple{4,Int64}=
                   (0,0,0,0))

    # BigWig output files
    mml1_path,mml2_path,nme1_path,nme2_path,uc_path = out_paths

    # Find chromosomes
    reader_fa = open(FASTA.Reader,fa,index=fa*".fai")
    chr_sizes = reader_fa.index.lengths
    chr_names = reader_fa.index.names
    close(reader_fa)

    # Loop over chromosomes
    for chr in chr_names

        # Get windows pertaining to current chromosome
        print_log("Processing chromosome $(chr) ...")
        haps_chr = read_gff_chr(gff,chr)
        chr_size = chr_sizes[findfirst(x->x==chr, chr_names)]

        # Process them
        out_pmap = pmap(hap -> proc_obs_hap(hap,chr,chr_size,bam1,bam2,gff,fa,out_paths,pe,g_max,
                        cov_ths,trim),haps_chr)

        # Add last to respective bedGraph file
        write_tobs([x[1] for x in out_pmap],chr,mml1_path)
        write_tobs([x[2] for x in out_pmap],chr,mml2_path)
        write_tobs([x[3] for x in out_pmap],chr,nme1_path)
        write_tobs([x[4] for x in out_pmap],chr,nme2_path)
        write_tobs([x[5] for x in out_pmap],chr,uc_path)

    end

    # Sort
    sort_bedgraphs([mml1_path,mml2_path,nme1_path,nme2_path,uc_path])

    # Return
    return nothing

end # end comp_tobs
"""
    `write_tnull(RECORDS,PATH)`

Function that writes null records in `RECORDS` into `PATH`.

# Examples
```julia-repl
julia> JuliASM.write_tnull(RECORDS,PATH)
```
"""
function write_tnull(recs::Vector{Tuple{Int64,Int64,Float64}},path::String)

    # Open output bedGraph file in append mode (no need to close it)
    open(path, "a") do f
        for i=1:length(recs)
            recs[i][1]!=0 || continue
            write(f,join(string.(collect(recs[i])),"\t"),"\n")
        end
    end

    # Return
    return nothing

end # end write_tnull
"""
    `get_kstar_table(GFF_PATH,CHR_NAMES,BLK_SIZE)`

Function that returns a table with maximum K for each N in GFF_PATH.

# Examples
```julia-repl
julia> JuliASM.get_kstar_table(GFF_PATH,CHR_NAMES,BLK_SIZE)
```
"""
function get_kstar_table(gff::String,chr_names::Vector{String},g_max::Int64)::Dict{Int64,Int64}

    # Loop over chromosomes
    table = Dict{Int64,Int64}()
    for chr in chr_names
        hom_win = read_gff_chr(gff,chr)
        # Loop over windows in chromosome
        for win in hom_win
            # Get window info
            win_st = GFF3.seqstart(win)
            win_end = GFF3.seqend(win)
            cpg_pos = get_cpg_pos(Dict(GFF3.attributes(win)))
            n = get_ns(cpg_pos[1],g_max,win_st,win_end)
            # Update table if necessary
            haskey(table,sum(n)) || (table[sum(n)]=length(n))
            length(n)>table[sum(n)] && (table[sum(n)]=length(n))
        end
    end

    # Return
    return table

end # end get_kstar_table
"""
    `subset_haps_cov(GFF,BAM,FA,PE,COV_THS,NTOT,TRIM,CHRmCHR_SIZE)`

Function that returns a set of haplotypes with N=Ntot and satisfying the coverage requirements.

# Examples
```julia-repl
julia> JuliASM.subset_haps_cov(gff,bam,fa,pe,cov_ths,ntot,trim,chr,chr_size)
```
"""
function subset_haps_cov(gff::String,bam::String,fa::String,pe::Bool,cov_ths::Int64,
                         n::Vector{Int64},trim::NTuple{4,Int64},chr::String,
                         chr_size::Int64)::Vector{GFF3.Record}

    # Ntot
    ntot = sum(n)

    # Load genomic features from a GFF3 file
    out_haps = Vector{GFF3.Record}()
    hom_haps = open(collect,GFF3.Reader,gff)

    # Filter based on N and chromosomes
    filter!(x -> (GFF3.score(x)>=ntot) && (GFF3.seqid(x) == chr),hom_haps)
    length(hom_haps)>0 || return out_haps

    # Loop over homozygous blocks
    for hap in hom_haps

        # All CpG sites in haplotype
        cpg_pos = get_cpg_pos(Dict(GFF3.attributes(hap)))[1]

        # Get ntot contiguous CpG sites at random
        ind = sample(1:(length(cpg_pos)-ntot+1))
        cpgs = cpg_pos[ind:(ind+ntot-1)]

        # Obtain observations overlapping with CpG sites
        xobs = read_bam(bam,chr,cpgs[1],cpgs[end],cpgs,chr_size,pe,trim)

        # Accept candidate haplotype if coverage is OK  in all subregions
            # if 4*cov_ths <= mean_cov(xobs) <= 400
        if all(2*cov_ths .<= mean_cov_sr(xobs,n) .<= 400)
            new_hap = "$(chr)\t.\t.\t$(cpgs[1])\t$(cpgs[end])\t$(ntot)\t.\t.\tCpGs=$(cpgs)"
            push!(out_haps,GFF3.Record(new_hap))
        end

        # Break if enough haps for chromosome
        length(out_haps)>500 && break

    end

    # Return
    return out_haps

end # end subset_haps_cov
"""
    `get_nvec_kstar(NTOT,KSTAR)`

Function that returns vector n with Ntot CpG sites divided into Kstar subregions.

# Examples
```julia-repl
julia> JuliASM.get_nvec_kstar(10,4)
4-element Array{Int64,1}:
 3
 2
 2
 3
```
"""
function get_nvec_kstar(ntot::Int64,kstar::Int64)::Vector{Int64}

    # Partition CpG sites into Kstar subregions
    n = div(ntot,kstar) * ones(Int64,kstar)
    if sum(n)<ntot
        cpg_ids = rand(Categorical(1.0/kstar*ones(kstar)),ntot-sum(n))
        n += [length(findall(x->x==i,cpg_ids)) for i=1:kstar]
    end

    # Return n
    return n

end # end get_nvec_kstar
"""
    `cov_obs_part(XOBS,N,COV_THS,COV_A,COV_B)`

Function that partitions observations into two sets of observations that satisfy a minimum
coverage within [COV_THS-COV_A,COV_THS+COV_B] for each subregion. If not partition is found that
satisfies it, then two empty vectors are returned.

# Examples
```julia-repl
julia> JuliASM.cov_obs_part(xobs,n,cov_ths,cov_a,cov_b)
```
"""
function cov_obs_part(xobs::Vector{Vector{Int64}},n::Vector{Int64},cov_ths::Int64,cov_a::Float64,
                      cov_b::Float64)::Tuple{Vector{Vector{Int64}},Vector{Vector{Int64}}}

    ## Part that takes care of minimum coverage. After this step an allele can have
    ## more than cov_ths+b. Here we just guarantee mean_cov(xobsi)≧cov_ths.

    # Divide into two sets of equal size while checking bot have at least cov_ths
    ct = 1
    fail = false
    xobs1 = Vector{Int64}()
    xobs2 = Vector{Int64}()
    while true
        # Partition
        ind = sample(1:length(xobs),div(length(xobs),2),replace=false)
        xobs1 = xobs[ind]
        xobs2 = xobs[setdiff(1:length(xobs),ind)]
        # If coverage okay break loop
        (mean_cov(xobs1)>=cov_ths) && (mean_cov(xobs2)>=cov_ths) && break
            # all(mean_cov_sr(xobs1,n).>=cov_ths) && all(mean_cov_sr(xobs2,n).>=cov_ths) &&  break
            # Check we haven't exceeded the maximum number of permutations
        if ct>20
            # We failed and break
            fail=true
            break
        else
            # Increase counter
            ct += 1
        end
    end

    # Return empty arrays if fail=true
    fail==false || return [],[]
        # println("First OK: $(mean_cov(xobs1)) & $(mean_cov(xobs2))")

    ## Part extra to limit number of observations used. The goal is to reduce the size of
    ## xobs1 & xobs2 until both coverages are within [cov_ths-cov_a,cov_ths+cov_b].

    # Check for [cov_ths-cov_a,cov_ths+cov_b] for xobs1
    less_ok = false
    while true
        # Check if we can keep decreasing number
        less_ok = mean_cov(xobs1[2:end])>=(cov_ths-cov_a) && mean_cov(xobs1)>=(cov_ths+cov_b)
            # less_ok = all(mean_cov_sr(xobs1[2:end],n).>=(cov_ths-cov_a)) &&
            # all(mean_cov_sr(xobs1,n).>=(cov_ths+cov_b))
        xobs1 = less_ok ? xobs1=xobs1[2:end] : xobs1
        less_ok || break
    end

    # Check if resulting observations meet b-bound as well
    fail = mean_cov(xobs1)<=(cov_ths+cov_b) ? false : true
        # fail = all(mean_cov_sr(xobs1,n).<=(cov_ths+cov_b)) ? false : true
    fail==false || return [],[]
        # println("Second OK: $(mean_cov(xobs1))")

    # Check for [cov_ths-cov_a,cov_ths+cov_b] for xobs2
    less_ok = false
    while true
        # Check if we can keep decreasing number
        less_ok = mean_cov(xobs2[2:end])>=(cov_ths-cov_a) && mean_cov(xobs2)>=(cov_ths+cov_b)
            # less_ok = all(mean_cov_sr(xobs2[2:end],n).>=(cov_ths-cov_a)) &&
            #     all(mean_cov_sr(xobs2,n).>=(cov_ths+cov_b))
        xobs2 = less_ok ? xobs2=xobs2[2:end] : xobs2
        less_ok || break
    end

    # Check if resulting observations meet b-bound as well
    fail = mean_cov(xobs2)<=(cov_ths+cov_b) ? false : true
        # fail = all(mean_cov_sr(xobs2,n).<=(cov_ths+cov_b)) ? false : true
    fail==false || return [],[]
        # println("Third OK: $(mean_cov(xobs2))")

    # Return partition
    return xobs1,xobs2

end # end cov_obs_part
"""
    `proc_null_hap()`

Function that computes a null dMML, dNME, and UC given a homozygous haplotype.

# Examples
```julia-repl
julia> proc_null_hap()
```
"""
function proc_null_hap(hap::GFF3.Record,ntot::Int64,bam::String,het_gff::String,hom_gff::String,
                       fa::String,kstar::Int64,out_paths::Vector{String},pe::Bool,g_max::Int64,
                       cov_ths::Int64,cov_a::Float64,cov_b::Float64,trim::NTuple{4,Int64},
                       mc_null::Int64,n_max::Int64,n_subset::Vector{Int64},chr_dic::Dict{String,
                       Int64})::Vector{Tuple{Float64,Float64,Float64}}

    # Empty output
    nan_out = (NaN,NaN,NaN)

    # Get homozygous window
    chr = GFF3.seqid(hap)
    chr_size = chr_dic[chr]

    # Sample ntot contiguous CpG sites and partition into Kstar
    cpg_pos = get_cpg_pos(Dict(GFF3.attributes(hap)))[1]
    n = get_nvec_kstar(ntot,kstar)

    # Check average coverage is within normal limits (watch for repetitive regions)
    xobs = read_bam(bam,chr,cpg_pos[1],cpg_pos[end],cpg_pos,chr_size,pe,trim)
    all(2*cov_ths .<= mean_cov_sr(xobs,n) .<= 400) || return nan_out
        # 2*cov_ths <= mean_cov(xobs) <= 400 || return nan_out

    # Obtain a number of null statistics with different permutations
    stats = Vector{Tuple{Float64,Float64,Float64}}()
    for i=1:10

        # Randomly partition observations (sample minimum coverage)
        xobs1,xobs2 = cov_obs_part(xobs,n,cov_ths,cov_a,cov_b)
        (length(xobs1)>0) && (length(xobs2)>0) || continue

        # Estimate each single-allele model and check if on boundary of parameter space
        θ1 = est_theta_sa(n,xobs1)
        check_boundary(θ1) && continue
        θ2 = est_theta_sa(n,xobs2)
        check_boundary(θ2) && continue

        # Estimate moments
        ∇1 = get_grad_logZ(n,θ1)
        ∇2 = get_grad_logZ(n,θ2)

        # Estimate output quantities
        nme1 = round(comp_nme_∇(n,θ1,∇1);digits=8)
        nme2 = round(comp_nme_∇(n,θ2,∇2);digits=8)
        uc = round(comp_uc(trues(ntot),trues(ntot),n,n,θ1,θ2,nme1,nme2);digits=8)

        # Append stats
        push!(stats,(round(abs(comp_mml_∇(n,∇1)-comp_mml_∇(n,∇2));digits=8),abs(nme1-nme2),uc))

    end

    # Return output
    return stats

end
"""
    `comp_tnull(BAM_PATH,GFF_PATH,FA_PATH,OUT_PATH)`

Function that computes null dMMLs, dNME, and UC from a BAM files at locations given by GFF that
contains the windows with not genetic variants and a FASTA file that contains the reference genome.

# Examples
```julia-repl
julia> comp_tnull(BAM_PATH,GFF_PATH,FA_PATH,OUT_PATH)
```
"""
function comp_tnull(bam::String,het_gff::String,hom_gff::String,fa::String,out_paths::Vector{String}
                    ;pe::Bool=true,g_max::Int64=300,cov_ths::Int64=5,cov_a::Float64=0.0,
                    cov_b::Float64=1.0,trim::NTuple{4,Int64}=(0,0,0,0),mc_null::Int64=5000,
                    n_max::Int64=20,n_subset::Vector{Int64}=collect(1:n_max))

    # BigWig output files
    dmml_path,dnme_path,uc_path = out_paths

    # Find relevant chromosomes and sizes
    reader_fa = open(FASTA.Reader,fa,index=fa*".fai")
    chr_dic = Dict(zip(reader_fa.index.names,reader_fa.index.lengths))
    chr_names = collect(keys(chr_dic))
    close(reader_fa)

    # Obtain Kstar for every N & maximum N to analyze
    print_log("Generating Kstar table ...")
    kstar_table = get_kstar_table(het_gff,collect(keys(chr_dic)),g_max)

    # Cap Nmax
    n_max = min(maximum(keys(kstar_table)),n_max)

    # Loop over Ns of interest
    for ntot in sort(collect(intersect(keys(kstar_table),n_subset)))

        # Skip if greater than n_max
        ntot>n_max && continue

        # Get kstar
        kstar = kstar_table[ntot]

        # Get n vector
        n = get_nvec_kstar(ntot,kstar)

        # Print current N
        out_dmml = Vector{Tuple{Int64,Int64,Float64}}()
        out_dnme = Vector{Tuple{Int64,Int64,Float64}}()
        out_uc = Vector{Tuple{Int64,Int64,Float64}}()

        # Store all haplotypes with enough coverage
        print_log("Scanning haplotypes ...")
        haps = vcat(pmap(chr -> subset_haps_cov(hom_gff,bam,fa,pe,cov_ths,n,trim,chr,chr_dic[chr]),
                         chr_names)...)

        # Keep in memory mc_null only for memory efficiency
        haps = length(haps)>mc_null ? haps[sample(1:length(haps),mc_null)] : haps

        # Sample windows with N=ntot
        print_log("Sampling a total of $(mc_null) null statistics ...")
        print_log("Candidate haplotypes with N=$(ntot) are $(length(haps)) ...")
        i = 1
        while length(out_dmml)<mc_null

            # Subsample
            haps_subset = length(haps)>1000 ? haps[sample(1:length(haps),1000)] : haps

            # Process them in parallel
            out_pmap = vcat(pmap(hap -> proc_null_hap(hap,ntot,bam,het_gff,hom_gff,fa,kstar,
                                 out_paths,pe,g_max,cov_ths,cov_a,cov_b,trim,mc_null,n_max,
                                 n_subset,chr_dic),haps_subset)...)

            # Keep only the ones with data
            out_pmap = out_pmap[map(stat->!any(isnan.(stat)),out_pmap)]

            # Save succesful runs
            if length(out_pmap)>0
                print_log("Saving $(length(out_pmap)) null stats ...")
                append!(out_dmml,map(stat->(ntot,kstar,stat[1]),out_pmap))
                append!(out_dnme,map(stat->(ntot,kstar,stat[2]),out_pmap))
                append!(out_uc,map(stat->(ntot,kstar,stat[3]),out_pmap))
            end

            # Break if run for too long
            if i>20
                print_log("Exceeded $(i-1) iterations in comp_tnull ...")
                break
            end

            # Increase counter
            i += 1

        end

        # Add last to respective bedGraph file
        write_tnull(out_dmml,dmml_path)
        write_tnull(out_dnme,dnme_path)
        write_tnull(out_uc,uc_path)

    end # end N for loop

    # Sort files
    sort_bedgraphs([dmml_path,dnme_path,uc_path])

    # Print message
    print_log("Done with null statistics ...")

    # Return nothing
    return nothing

end # end comp_tnull
"""
    `comp_pvals_stat(TOBS_PATH,TNULL_PATH,N_MAX)`

Function that computes adjusted p-values for each haplotype in TOBS_PATH from the empirical
distribution of the pertaining statistic in TNULL_PATH. The output file has the genomic
coordinates, the value of the statistic, and the adjusted p-value.

# Examples
```julia-repl
julia> JuliASM.comp_pvals_stat(tobs_path,tnull_path)
```
"""
function comp_pvals_stat(tobs_path::Vector{String},tnull_path::String,p_path::String,n_max::Int64)

    # Get null stats
    tnull = readdlm(tnull_path,'\t',Any)

    # Consider case dMML/dNME VS UC
    if length(tobs_path)>1
        # dMML/dNME. NOTE: both files need to have the same order of haplotypes!
        tobs = readdlm(tobs_path[1],'\t',Any)
        tobs[:,4] = abs.(tobs[:,4]-readdlm(tobs_path[2],'\t',Any)[:,4])
    else
        # UC
        tobs = readdlm(tobs_path[1],'\t',Any)
    end

    # Initialize p-value matrix
    pvals = tobs[:,1:5]
    pvals[:,5] .= "NO NULL DATA"

    # Compute p-values
    for i=1:size(tobs)[1]
        # Check we can compute p-value
        (tobs[i,5]<=n_max && sum(tnull[:,1].==tobs[i,5])>0) || continue
        # Compute empirical p-value
        pvals[i,5] = sum(tnull[tnull[:,1].==tobs[i,5],3].>=tobs[i,4]) / sum(tnull[:,1].==tobs[i,5])
    end

    # Multiple hypothesis testing correction. NOTE: should we apply BH on each N independently?
    ind = pvals[:,5].!="NO NULL DATA"
    pvals[ind,5] = MultipleTesting.adjust(convert(Vector{Float64},pvals[ind,5]),BenjaminiHochberg())

    # Write to output
    open(p_path,"w") do io
        writedlm(io,pvals,'\t')
    end

    # Return
    return nothing

end # end comp_pvals_stat
"""
    `comp_pvals(TOBS_PATH,TNULL_PATH,PVALS_PATH,N_MAX)`

Function that computes p-values from empirical null distribution.

# Examples
```julia-repl
julia> comp_pvals(tobs_path,tnull_path,p_path,n_max)
```
"""
function comp_pvals(tobs_path::Vector{String},tnull_path::Vector{String},p_path::Vector{String},
                    n_max::Int64)

    # Compute three sets of pvalues
    print_log("Computing p-values dMML ...")
    comp_pvals_stat(tobs_path[1:2],tnull_path[1],p_path[1],n_max)
    print_log("Computing p-values dNME ...")
    comp_pvals_stat(tobs_path[3:4],tnull_path[2],p_path[2],n_max)
    print_log("Computing p-values UC ...")
    comp_pvals_stat([tobs_path[5]],tnull_path[3],p_path[3],n_max)
    print_log("Done with p-values ...")

    # Return
    return nothing

end # end comp_pvals
"""
    `run_analysis(BAM1_PATH,BAM2_PATH,BAMU_PATH,VCF_PATH,FA_PATH,OUT_PATH)`

Function that estimates MML1/2, NME1/2, and NMI on a pair of BAM files (allele 1, allele 2) given a
(phased) VCF file that contains the heterozygous SNPs and a FASTA file that contains the reference
genome. An empirical null distribution for each quantity is estimated from a set of representative
homozygous regions from unassigned BAM records in BAMU_PATH. A p-value is computed for each
haplotype in the VCF file.

# Examples
```julia-repl
julia> run_analysis(BAM1_PATH,BAM2_PATH,BAMU_PATH,VCF_PATH,FA_PATH,OUT_PATH)
```
"""
function run_analysis(bam1::String,bam2::String,bamu::String,vcf::String,fa::String,outdir::String;
                      pe::Bool=true,g_max::Int64=300,win_exp::Int64=100,cov_ths::Int64=5,
                      cov_a::Float64=0.0,cov_b::Float64=1.0,trim::NTuple{4,Int64}=(0,0,0,0),
                      mc_null::Int64=5000,n_max::Int64=25,n_subset::Vector{Int64}=collect(1:n_max))

    # Print initialization of juliASM
    print_log("Starting JuliASM analysis ...")

    # Check index files exist
    if !(isfile.(bam1*".bai",bam1*".bai",fa*".fai"))
        print_log("Index files for BAM or FASTA missing. Exiting julia ...")
        exit(1)
    end

    # Create output folder if it doesn't exist
    isdir(outdir) || mkdir(outdir)

    # BigWig output files
    prefix_sample = String(split(basename(bam1),".")[1])
    tobs_path,tnull_path,p_path = get_outpaths(outdir,prefix_sample)

    # Check for existance of at least an output files
    if all(isfile.(vcat(tobs_path,tnull_path)))
        print_log("At least an output file already exists. Exiting julia ...")
        exit(1)
    end

    # Create gff file with heterozygous loci
    print_log("Reading in VCF & FASTA files ...")
    het_gff = outdir * split(basename(vcf),".")[1] * "_het.juliasm.gff"
    hom_gff = outdir * split(basename(vcf),".")[1] * "_hom.juliasm.gff"
    if !(isfile(het_gff))
        print_log("Generating JuliASM GFF files ...")
        gen_gffs([het_gff,hom_gff],fa,vcf,win_exp,n_max)
    end

    # Compute observed statistics from heterozygous loci
    print_log("Computing observed statistics in heterozygous loci...")
    comp_tobs(bam1,bam2,het_gff,fa,tobs_path;pe=pe,g_max=g_max,cov_ths=cov_ths,trim=trim)

    # Compute null statistics from homozygous loci
    print_log("Generating null statistics in homozygous loci ...")
    comp_tnull(bamu,het_gff,hom_gff,fa,tnull_path;pe=pe,g_max=g_max,cov_ths=cov_ths,
               cov_a=cov_a,cov_b=cov_b,trim=trim,mc_null=mc_null,n_max=n_max,n_subset=n_subset)

    # Compute null statistics from heterozygous loci
    print_log("Computing p-values in heterozygous loci ...")
    comp_pvals(tobs_path,tnull_path,p_path,n_max)

    # Print number of features interrogated and finished message
    print_log("Done.")

    # Return
    return nothing

end # end run_analysis
