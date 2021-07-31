##################################################
# ASVs

rule dada2_download_taxa:
    output:
        os.path.join(DBS_DIR, "dada2", os.path.basename(config["dada2"]["asv"]["tax"]))
    params:
        url=config["dada2"]["asv"]["tax"]
    log:
        os.path.join(DBS_DIR, "dada2", os.path.basename(config["dada2"]["asv"]["tax"]) + ".log")
    message:
        "DADA2: download taxa data"
    shell:
        "(date && wget -O {output} {params.url} && date) &> {log}"

rule dada2_asv_errormod:
    input:
        r1=expand("fastq/cutadapt_adapters/{sid}_R1.fastq.gz", sid=SAMPLES),
        r2=expand("fastq/cutadapt_adapters/{sid}_R2.fastq.gz", sid=SAMPLES)
    output:
        pdf="dada2/errormod.pdf",
        rds="dada2/errormod.rds"
    log:
        "dada2/errormod.log"
    threads:
        config["dada2"]["asv"]["threads"]
    params:
        utils=os.path.join(SRC_DIR, "utils.R")
    conda:
        os.path.join(ENV_DIR, "r.yaml")
    message:
        "DADA2: error model"
    script:
        os.path.join(SRC_DIR, "dada2_asv_errormod.R")

rule dada2_asv:
    input:
        r1="fastq/cutadapt_adapters/{sid}_R1.fastq.gz",
        r2="fastq/cutadapt_adapters/{sid}_R2.fastq.gz",
        errormod="dada2/errormod.rds"
    output:
        rds=temp("dada2/{sid}.rds")
    log:
        "dada2/{sid}.log"
    threads:
        config["dada2"]["asv"]["threads_single"]
    params:
        utils=os.path.join(SRC_DIR, "utils.R")
    conda:
        os.path.join(ENV_DIR, "r.yaml")
    message:
        "DADA2: ASVs: {input}"
    script:
        os.path.join(SRC_DIR, "dada2_asv.R")

rule dada2_merge:
    input:
        asv=expand("dada2/{sid}.rds", sid=SAMPLES),
        tax=os.path.join(DBS_DIR, "dada2", os.path.basename(config["dada2"]["asv"]["tax"]))
    output:
        counts="dada2/ASV.counts.tsv",
        fasta="dada2/ASV.fasta",
        tax="dada2/ASV.tax.tsv",
        stats="dada2/ASV.stats.tsv",
        rds="dada2/ASV.rds"
    log:
        "dada2/ASV.log"
    threads:
        config["dada2"]["asv"]["threads"]
    params:
        utils=os.path.join(SRC_DIR, "utils.R")
    conda:
        os.path.join(ENV_DIR, "r.yaml")
    message:
        "DADA2: merge"
    script:
        os.path.join(SRC_DIR, "dada2_asv_merge.R")    

##################################################
# Stats

rule fastq_asv_stats:
    input:
        raw="multiqc/fastqc/raw/multiqc_data/multiqc_fastqc.txt",
        trim1="multiqc/fastqc/cutadapt_primers/multiqc_data/multiqc_fastqc.txt",
        trim2="multiqc/fastqc/cutadapt_adapters/multiqc_data/multiqc_fastqc.txt",
        asv="dada2/ASV.stats.tsv"
    output:
        tsv="dada2/stats.tsv",
        pdf="dada2/stats.pdf"
    log:
        "dada2/stats.log"
    params:
        utils=os.path.join(SRC_DIR, "utils.R")
    conda:
        os.path.join(ENV_DIR, "r.yaml")
    message:
        "FastQ/ASV stats"
    script:
        os.path.join(SRC_DIR, "fastq_asv_stats.R")


##################################################
# Krona

# Krona text input: https://github.com/marbl/Krona/wiki/Importing-text-and-XML-data
rule krona_text:
    input:
        counts="dada2/ASV.counts.tsv",
        tax="dada2/ASV.tax.tsv"
    output:
        temp(expand("dada2/ASV.krona.{sid}.txt", sid=SAMPLES))
    threads:
        1
    message:
        "KronaTools: TEXT input from {input}"
    run:
        import re
        from pandas import read_csv
        # read in data
        counts = read_csv(input.counts, sep='\t', header=0, index_col=0) # sample x asv
        tax    = read_csv(input.tax,    sep='\t', header=0, index_col=0) # asv x taxonomy
        # create output files
        for sid in counts.index:
            sid_re = re.compile(".*\.krona\.%s\.txt$" % sid) # pattern should match the output file pattern
            with open(list(filter(sid_re.match, output))[0], "w") as ofile: # matching output file
                for asv in counts.columns:
                    if counts.loc[sid, asv] > 0:
                        assert asv in tax.index
                        c = counts.loc[sid, asv] # count
                        t = "\t".join(tax.loc[asv].fillna(value="NA")) # taxonomy
                        t = re.sub("^NA(\tNA)+$", "Unknown", t) # replace completely unknown taxonomy
                        t = re.sub("(\tNA)+$", "", t) # remove trailing NAs
                        ofile.write("%d\t%s\n" % (c, t))

rule krona_textimport:
    input:
        expand("dada2/ASV.krona.{sid}.txt", sid=SAMPLES)
    output:
        "dada2/ASV.krona.html"
    threads:
        1
    conda:
        os.path.join(ENV_DIR, "krona.yaml")
    message:
        "KronaTools: plot from TEXT input: {input}"
    shell:
        "ktImportText {input} -o {output} && sed -i 's/ASV\.krona\.//g' {output}"

##################################################
# Phyloseq

rule phyloseq:
    input:
        counts="dada2/ASV.counts.tsv",
        tax="dada2/ASV.tax.tsv",
        meta=config["samples_meta"]
    output:
        pdf="dada2/ASV.phyloseq.pdf"
    log:
        "dada2/ASV.phyloseq.log"
    params:
        utils=os.path.join(SRC_DIR, "utils.R"),
        ord_method="DCA",
        ord_dist="bray",
    conda:
        os.path.join(ENV_DIR, "r.yaml")
    message:
        "Phyloseq"
    script:
        os.path.join(SRC_DIR, "phyloseq.R") 