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
julia> try_olaps(reader,chr,win)
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
julia> read_bam(BAM_PATH,"chr1",30,80,[40,60],false,(0,0,0,0))
```
"""
function read_bam(bam::String,chr::String,f_st::Int64,f_end::Int64,cpg_pos::Vector{Int64},
                  chr_size::Int64,pe::Bool,trim::NTuple{4,Int64})::Array{Vector{Int64},1}

    # Number of CpG sites is determined by that in the region
    N = length(cpg_pos)

    # Expand window in PE case to include pair even if it's outside the window
    exp_win = pe ? [max(1,f_st-75),min(chr_size,f_end+75)] : [f_st,f_end]

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
julia> write_gff!(gff, gff_records)
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

end # end write_gff
"""
    `next_record(READER,CHR_NAMES,RECORD)`

Recursive function that returns next record in stream if in chr_names, and false if the stream
reached its end. This function is only used in `gen_gffs()`.

# Examples
```julia-repl
julia> next_record(reader, chr_names, record)
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
julia> is_het_cpg!(var,seq,h1,h2)
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

end # end is_het_cpg!
"""
    `get_records_ps!(VCF_READER,SEQ,RECORD,CURR_PS,WIN_END,H1,H2)`

Recursive function that returns the end of current PS phased SNPs along the next SNP. This
function is only used in `gen_gffs()`.

# Examples
```julia-repl
julia> get_records_ps!(reader,seq,record,curr_ps,wEnd,h1,h2)
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
            record,wEnd=get_records_ps!(reader,seq,record,curr_ps,wEnd,h1,h2)
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
julia> expand_win(50,80,10,50,2000)
??
```
"""
function expand_win(l_snp::Int64,r_snp::Int64,win_exp::Int64,blk_size::Int64,
                    chr_size::Int64)::Vector{Int64}

    # Expand left and right to include nearby CpG sites
    l_snp -= win_exp
    r_snp += win_exp

    # Add remaining to be multiple of blk_size
    rmdr = (r_snp-l_snp)%blk_size==0 ? 0 : div(blk_size-(r_snp-l_snp)%blk_size,2)

    # Cap at 1 and end of chromosome
    return [max(1,l_snp-rmdr),min(chr_size,r_snp+rmdr)]

end # end expand_win
"""
    `find_num_models(CPG_POS,NMAX,NMOD)`

Function that returns the number of models required so as to avoid joint model with N>NMAX.

# Examples
```julia-repl
julia> find_num_models([50,55,60,65],2,1)
2
```
"""
function find_num_models(cpg_pos::Vector{Int64},n_max::Int64,n_mod::Int64)::Int64

    # Check if resulting models are acceptable
    return length(cpg_pos)/n_mod>n_max ? find_num_models(cpg_pos,n_max,n_mod+1) : n_mod

end # end find_num_models
"""
    `gen_gffs([HET_GFF_PATH,HOM_GFF_PATH],FA_PATH,VCF_PATH,WIN_EXP,BLOCK_SIZE,NMAX)`

Function that creates a GFF file containing the genomic regions to be analyzed by `JuliASM`. If a
set of SNPs is phased, this function will create a single window for all. In case, the set of SNPs
is not phased, then an individual window will be created for each SNP. The phasing of the VCF
records should be specified by means of the standard `PS` tag. In addition, a GFF file containing
the homozygous portion of the DNA is also generated together with the position of the CpG sites
in it.

# Examples
```julia-repl
julia> gen_gffs([HET_GFF_PATH,HOM_GFF_PATH],FA_PATH,VCF_PATH,WIN_EXP,BLOCK_SIZE,NMAX)
```
"""
function gen_gffs(gff::Vector{String},fa::String,vcf::String,win_exp::Int64,blk_size::Int64,
                  n_max::Int64)

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
        record = next_record(reader_vcf, chr_names, record)

        # Check we have loaded the right chromosome from FASTA
        if !(curr_chr == VCF.chrom(record))
            curr_chr = VCF.chrom(record)
            reader_fa = open(FASTA.Reader,fa,index=fa*".fai")
            fa_record = reader_fa[curr_chr]
            chr_size = chr_sizes[findfirst(x->x==curr_chr, chr_names)]
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

        # Get remainder from block size and distribute it evenly
        win = expand_win(wSt,wEnd,win_exp,blk_size,chr_size)

        # Obtain DNA sequence from reference and homozygous CpG loci from it
        wSeq = convert(String,FASTA.sequence(fa_record,win[1]:win[2]))
        cpg_pos = map(x->getfield(x,:offset),eachmatch(r"CG",wSeq)).+win[1].-1

        # Store heterozygous haplotype & corresponding CpG sites if homozygous CpG site/s
        if length(cpg_pos)>0

            # Obtain required number of models
            n_mod = find_num_models(cpg_pos,n_max,1)
            mod_size = div(length(cpg_pos),n_mod)
            mod_rem = length(cpg_pos)-n_mod*mod_size

            # Print each set separately
            for i=1:n_mod
                # Get the mod_size homozygous CpG sites
                cpg_pos_mod = cpg_pos[(mod_size*(i-1)+1):min((mod_size*i),length(cpg_pos))]
                i==n_mod && mod_rem>0 && append!(cpg_pos_mod,cpg_pos[(end-mod_rem+1):end])
                win_mod = [cpg_pos_mod[1],cpg_pos_mod[end]]
                # Add heterozygous part
                het1_mod = [in(het1[i],cpg_pos_mod[1]:cpg_pos_mod[end]) for i=1:length(het1)]
                het2_mod = [in(het2[i],cpg_pos_mod[1]:cpg_pos_mod[end]) for i=1:length(het2)]
                length(het1_mod)>0 && i==1 && minimum(het1)<cpg_pos[1] && (het1_mod[1]=true)
                length(het2_mod)>0 && i==1 && minimum(het2)<cpg_pos[1] && (het2_mod[1]=true)
                length(het1_mod)>0 && i==n_mod && maximum(het1)>cpg_pos[end] && (het1_mod[end]=true)
                length(het2_mod)>0 && i==n_mod && minimum(het2)>cpg_pos[end] && (het2_mod[end]=true)
                i==1 && length(het1)>0 && minimum(het1)<cpg_pos[1] && (win_mod[1]=minimum(het1))
                i==1 && length(het2)>0 && minimum(het2)<cpg_pos[1] && (win_mod[1]=minimum(het2))
                i==n_mod && length(het1)>0 && maximum(het1)>cpg_pos[end] &&
                    (win_mod[end]=maximum(het1))
                i==n_mod && length(het2)>0 && maximum(het2)>cpg_pos[end] &&
                    (win_mod[end]=maximum(het2))
                # Construct output string
                out_str = "$(curr_chr)\t.\t$(curr_ps)\t$(win_mod[1])\t$(win_mod[2])\t.\t.\t.\t"
                out_str *= "N=$(length(cpg_pos_mod));CpGs=$(cpg_pos_mod)"
                length(het1_mod)>0 && sum(het1_mod)>0 && (out_str*=";hetCpGg1=$(het1[het1_mod])")
                length(het2_mod)>0 && sum(het2_mod)>0 && (out_str*=";hetCpGg2=$(het2[het2_mod])")
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

end # end gen_gffs
"""
    `read_gff_chr(GFF_PATH,CHR)`

Function that reads in a GFF3 file in `GFF_PATH` and returns an array of GenomicFeatures Record
objects contained in chromosome `CHR`. The format is the standard defined by ENSEMBL in
https://useast.ensembl.org/info/website/upload/gff3.html.

# Examples
```julia-repl
julia> read_gff_chr(GFF_PATH,"chr1")
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
julia> write_tobs(RECORDS,CHR,PATH)
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

end # end write_tobs
"""
    `get_cpg_pos(FEAT_ATTS)`

Function that takes the attributes in the GFF file and returns the homozygous and heterozygous CpG
sites in it.

# Examples
```julia-repl
julia> get_cpg_pos(feat_atts)
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
    `get_ns(CPG_POS,BLOCK_SIZE,FEAT_ST)`

Functions that returns `[N1,...,NK]` given the position of the CpG sites and a block size
`BLOCK_SIZE` into which the region is partitioned.

# Examples
```julia-repl
julia> get_ns([100,250,300,350],200,90)
2-element Vector{Int64}:
 2
 2
```
"""
function get_ns(cpg_pos::Vector{Int64},blk_size::Int64,f_st::Int64)::Vector{Int64}

    # Check all cpg_pos are not upstream of f_st (provisional check)
    if minimum(cpg_pos)<f_st
        println(stderr, "[$(now())]: CpG sites in $cpg_pos upstream of f_st at $f_st.")
        println(stderr, "[$(now())]: Exiting julia ...")
        exit(1)
    end

    # Partition array
    y = Int64[]
    for p in cpg_pos
        j = 0
        while true
            ((p-f_st-j*blk_size)>=0) && ((p-f_st-(j+1)*blk_size)<=0) && break
            j += 1
        end
        push!(y,j)
    end

    # Return array of arrays
    return [sum(y.==i) for i in unique(y)]

end # end get_ns
"""
    `mean_cov(XOBS)`

Function returns the average coverage per CpG given some observations `XOBS`.

# Examples
```julia-repl
julia> xobs=[[1,-1] for i=1:10]; append!(xobs,[[1,0] for i=1:10]);
julia> mean_cov(xobs)
15.0
```
"""
function mean_cov(xobs::Array{Vector{Int64},1})::Float64

    # Return 0 if no observations
    return length(xobs)>0 ? norm(hcat(xobs...),1)/length(xobs[1]) : 0.0

end # end mean_cov
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

    # Return path for Tobs and Tnull
    return tobs_path,tnull_path

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

end # end print_log
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
                   pe::Bool=true,blk_size::Int64=300,cov_ths::Int64=5,trim::NTuple{4,Int64}=
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
        print_log("Processing 🧬  $(chr) ...")
        features_chr = read_gff_chr(gff,chr)
        chr_size = chr_sizes[findfirst(x->x==chr, chr_names)]

        # bedGraph records
        out_sa = SharedArray{Tuple{Int64,Int64,Float64,Int64,Int64},2}(length(features_chr),8)

        # Loop over windows in chromosome
        @sync @distributed for i=1:length(features_chr)

            # Get window of interest
            feat = features_chr[i]
            f_st = GFF3.seqstart(feat)
            f_end = GFF3.seqend(feat)

            # Get homozygous & heterozygous CpG sites (1:hom, 2:hap1, 3: hap2)
            cpg_pos = get_cpg_pos(Dict(GFF3.attributes(feat)))
            n = get_ns(cpg_pos[1],blk_size,f_st)
            n1 = get_ns(cpg_pos[2],blk_size,f_st)
            n2 = get_ns(cpg_pos[3],blk_size,f_st)

            # Get vectors from BAM1/2 overlapping feature & compute average coverage
            xobs1 = read_bam(bam1,chr,f_st,f_end,cpg_pos[2],chr_size,pe,trim)
            mean_cov1 = mean_cov(xobs1)
            mean_cov1>=cov_ths || continue
            xobs2 = read_bam(bam2,chr,f_st,f_end,cpg_pos[3],chr_size,pe,trim)
            mean_cov2 = mean_cov(xobs2)
            mean_cov2>=cov_ths || continue

            # Estimate each single-allele model and check if on boundary of parameter space
            theta1 = est_theta(n1,xobs1)
            check_boundary(theta1) && continue
            theta2 = est_theta(n2,xobs2)
            check_boundary(theta2) && continue

            # Estimate first and second moments
            ex1 = comp_ex(n1,theta1[1:(end-1)],theta1[end])
            ex2 = comp_ex(n2,theta2[1:(end-1)],theta2[end])
            exx1 = comp_exx(n1,theta1[1:(end-1)],theta1[end])
            exx2 = comp_exx(n2,theta2[1:(end-1)],theta2[end])

            # Add results to output array
            out_sa[i,1] = (f_st,f_end,comp_mml(ex1),sum(n1),length(n1))
            out_sa[i,2] = (f_st,f_end,comp_mml(ex2),sum(n2),length(n2))
            out_sa[i,3] = (f_st,f_end,comp_nme(trues(sum(n1)),n1,theta1[1:(end-1)],theta1[end],
                ex1,exx1),sum(n1),length(n1))
            out_sa[i,4] = (f_st,f_end,comp_nme(trues(sum(n2)),n2,theta2[1:(end-1)],theta2[end],
                ex2,exx2),sum(n2),length(n2))

            # Add UC to output
            z1 = BitArray([p in cpg_pos[1] ? true : false for p in cpg_pos[2]])
            z2 = BitArray([p in cpg_pos[1] ? true : false for p in cpg_pos[3]])
            nme1 = comp_nme(z1,n1,theta1[1:(end-1)],theta1[end],ex1,exx1)
            nme2 = comp_nme(z2,n2,theta2[1:(end-1)],theta2[end],ex2,exx2)
            out_sa[i,5] = (f_st,f_end,comp_uc(z1,z2,n1,n2,theta1,theta2,nme1,nme2),sum(n),length(n))

        end

        # Add last to respective bedGraph file
        write_tobs(out_sa[:,1],chr,mml1_path)
        write_tobs(out_sa[:,2],chr,mml2_path)
        write_tobs(out_sa[:,3],chr,nme1_path)
        write_tobs(out_sa[:,4],chr,nme2_path)
        write_tobs(out_sa[:,5],chr,uc_path)

    end

    # Sort
    sort_bedgraphs([mml1_path,mml2_path,nme1_path,nme2_path,uc_path])

end # end comp_tobs
"""
    `write_tnull(RECORDS,PATH)`

Function that writes null records in `RECORDS` into `PATH`.

# Examples
```julia-repl
julia> write_tnull(RECORDS,PATH)
```
"""
function write_tnull(recs::Vector{Tuple{Int8,Int8,Float64}},path::String)

    # Open output bedGraph file in append mode (no need to close it)
    open(path, "a") do f
        for i=1:length(recs)
            recs[i][1]!=0 || continue
            write(f,join(string.(collect(recs[i])),"\t"),"\n")
        end
    end

end # end write_tnull
"""
    `get_kstar_table(GFF_PATH,CHR_NAMES,BLK_SIZE)`

Function that returns a table with maximum K for each N in GFF_PATH.

# Examples
```julia-repl
julia> get_kstar_table(GFF_PATH,CHR_NAMES,BLK_SIZE)
```
"""
function get_kstar_table(gff::String,chr_names::Vector{String},blk_size::Int64)::Dict{Int64,Int64}

    # Loop over chromosomes
    table = Dict{Int64,Int64}()
    for chr in chr_names
        hom_win = read_gff_chr(gff,chr)
        # Loop over windows in chromosome
        for win in hom_win
            # Get window info
            win_st = GFF3.seqstart(win)
            cpg_pos = get_cpg_pos(Dict(GFF3.attributes(win)))
            n = get_ns(cpg_pos[1],blk_size,win_st)
            # Update table if necessary
            haskey(table,sum(n)) || (table[sum(n)]=length(n))
            length(n)>table[sum(n)] && (table[sum(n)]=length(n))
        end
    end

    return table

end # end get_kstar_table
"""
    `sample_win_n(GFF_PATH,CHR_NAMES,NTOT,MC_NULL)`

Function that returns a set of size MC_NULL of homozygous genomic windows with at least NTOT CpG
sites.

# Examples
```julia-repl
julia> sample_win_n(GFF_PATH,CHR_NAMES,NTOT,MC_NULL)
```
"""
function sample_win_ntot(gff::String,chr_names::Vector{String},ntot::Int64,
                         mc_null::Int64)::Vector{GFF3.Record}

    # Load genomic features from a GFF3 file
    features = open(collect,GFF3.Reader,gff)

    # Filter based on N and chromosomes
    filter!(x -> (GFF3.score(x)>=ntot) && (GFF3.seqid(x) in chr_names),features)

    # Return
    return features[sample(1:length(features),mc_null;replace=true)]

end # end sample_win_n
"""
    `sample_ntot_cpgs(NTOT,CPG_POS)`

Function that samples NTOT contiguous CpG sites from CPG_POS.

# Examples
```julia-repl
julia> sample_ntot_cpgs(4,[10,15,20,25,30,35,40])
4-element Array{Int64,1}:
 10
 15
 20
 25
```
"""
function sample_ntot_cpgs(ntot::Int64,cpg_pos::Vector{Int64})::Vector{Int64}

    # Get start index
    st_ind = length(cpg_pos)>ntot ? sample(1:(length(cpg_pos)-(ntot-1))) : 1

    # Return n
    return cpg_pos[st_ind:(st_ind+ntot-1)]

end # end sample_ntot_cpgs
"""
    `get_nvec_kstar(NTOT,KSTAR)`

Function that returns vector n with Ntot CpG sites divided into Kstar subregions.

# Examples
```julia-repl
julia> get_nvec_kstar(10,4)
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
    `comp_tnull(BAM_PATH,GFF_PATH,FA_PATH,OUT_PATH)`

Function that computes null dMMLs, dNME, and UC from a BAM files at locations given by GFF that
contains the windows with not genetic variants and a FASTA file that contains the reference genome.

# Examples
```julia-repl
julia> comp_tnull(BAM_PATH,GFF_PATH,FA_PATH,OUT_PATH)
```
"""
function comp_tnull(bam::String,het_gff::String,hom_gff::String,fa::String,out_paths::Vector{String}
                    ;pe::Bool=true,blk_size::Int64=300,cov_ths::Int64=5,trim::NTuple{4,Int64}=
                    (0,0,0,0),mc_null::Int64=5000,n_max::Int64=20)

    # BigWig output files
    dmml_path,dnme_path,uc_path = out_paths

    # Find relevant chromosomes and sizes
    reader_fa = open(FASTA.Reader,fa,index=fa*".fai")
    chr_dic = Dict(zip(reader_fa.index.names,reader_fa.index.lengths))
    close(reader_fa)

    # Obtain Kstar for every N & maximum N to analyze
    print_log("Generating Kstar table ...")
    kstar_table = get_kstar_table(het_gff,collect(keys(chr_dic)),blk_size)

    # Cap Nmax
    n_max = min(maximum(keys(kstar_table)),n_max)

    # Loop over Ns of interes
    for ntot=1:n_max

        # Print current N
        print_log("Generating null statistics for N=$(ntot) ...")

        # Get kstar
        kstar = kstar_table[ntot]

        # Sample windows with N=ntot
        print_log("Sampling $(mc_null) windows ...")
        win_ntot = sample_win_ntot(hom_gff,collect(keys(chr_dic)),ntot,mc_null)

        # bedGraph records
        print_log("Initializing output array ...")
        out_sa = SharedArray{Tuple{Int8,Int8,Float64},2}(mc_null,3)

        # Produce a set of null statistics for each index
        print_log("Computing null statistis ...")
        @sync @distributed for i=1:length(win_ntot)

            # Get homozygous window
            feat = win_ntot[i]
            chr = GFF3.seqid(feat)
            chr_size = chr_dic[chr]

            # Sample ntot contiguous CpG sites and partition into Kstar
            cpg_pos = get_cpg_pos(Dict(GFF3.attributes(feat)))[1]
            cpg_pos = sample_ntot_cpgs(ntot,cpg_pos)
            n = get_nvec_kstar(ntot,kstar)

            # Check average coverage
            xobs = read_bam(bam,chr,cpg_pos[1],cpg_pos[end],cpg_pos,chr_size,pe,trim)
            mean_cov(xobs)>=2*cov_ths || continue

            # Randomly partition observations (sample minimum coverage)
            part = sample(1:(2*cov_ths),cov_ths,replace=false)
            xobs1 = xobs[part]
            xobs2 = xobs[setdiff(1:length(xobs),part)]

            # Estimate each single-allele model and check if on boundary of parameter space
            theta1 = est_theta(n,xobs1)
            check_boundary(theta1) && continue
            theta2 = est_theta(n,xobs2)
            check_boundary(theta2) && continue

            # Estimate moments
            ex1 = comp_ex(n,theta1[1:(end-1)],theta1[end])
            ex2 = comp_ex(n,theta2[1:(end-1)],theta2[end])
            exx1 = comp_exx(n,theta1[1:(end-1)],theta1[end])
            exx2 = comp_exx(n,theta2[1:(end-1)],theta2[end])

            # Estimate output quantities
            mml1 = comp_mml(ex1)
            mml2 = comp_mml(ex2)
            nme1 = comp_nme(trues(ntot),n,theta1[1:(end-1)],theta1[end],ex1,exx1)
            nme2 = comp_nme(trues(ntot),n,theta2[1:(end-1)],theta2[end],ex2,exx2)
            uc = comp_uc(trues(ntot),trues(ntot),n,n,theta1,theta2,nme1,nme2)

            # Store output
            out_sa[i,1] = (ntot,kstar,abs(round(mml1-mml2;digits=8)))
            out_sa[i,2] = (ntot,kstar,abs(round(nme1-nme2;digits=8)))
            out_sa[i,3] = (ntot,kstar,uc)

        end

        # Add last to respective bedGraph file
        write_tnull(out_sa[:,1],dmml_path)
        write_tnull(out_sa[:,2],dnme_path)
        write_tnull(out_sa[:,3],uc_path)

    end

    # Sort files
    sort_bedgraphs([dmml_path,dnme_path,uc_path])

end # end comp_tnull
"""
    `comp_pvals(TOBS_PATH,TNULL_PATH)`

Function that computes p-values from empirical null distribution.

# Examples
```julia-repl
julia> comp_pvals(tobs_path,tnull_path)
??
```
"""
function comp_pvals(tobs_path::Vector{String},tnull_path::Vector{String},n_max::Int64)

    # Loop over N
    for n=1:n_max

        # Store empirical in table (value,count)

        # Get p-val for every haplotype with N=n (parallelization)

    end

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
                      pe::Bool=true,blk_size::Int64=300,win_exp::Int64=100,cov_ths::Int64=5,
                      trim::NTuple{4,Int64}=(0,0,0,0),mc_null::Int64=5000,n_max::Int64=25)

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
    tobs_path,tnull_path = get_outpaths(outdir,prefix_sample)

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
        gen_gffs([het_gff,hom_gff],fa,vcf,win_exp,blk_size,n_max)
    end

    # Compute observed statistics from heterozygous loci
    print_log("Computing observed statistics in heterozygous loci...")
    comp_tobs(bam1,bam2,het_gff,fa,tobs_path;pe=pe,blk_size=blk_size,cov_ths=cov_ths,trim=trim)

    # Compute null statistics from homozygous loci
    print_log("Generating null statistics in homozygous loci ...")
    comp_tnull(bamu,het_gff,hom_gff,fa,tnull_path;pe=pe,blk_size=blk_size,cov_ths=cov_ths,
        trim=trim,mc_null=mc_null,n_max=n_max)

    # Compute null statistics from heterozygous loci
    print_log("Computing p-values in heterozygous loci ...")
    comp_pvals(tobs_path,tnull_path,n_max)

    # Print number of features interrogated and finished message
    print_log("Done.")

end # end run_analysis
