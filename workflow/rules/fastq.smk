##################################################
# FASTQ files

localrules: link_fastq

rule link_fastq:
    input:
        r1=lambda wildcards: config["samples"][wildcards.sid]["R1"],
        r2=lambda wildcards: config["samples"][wildcards.sid]["R2"]
    output:
        r1="fastq/raw/{sid}_R1.fastq.gz",
        r2="fastq/raw/{sid}_R2.fastq.gz"
    wildcard_constraints:
        sid="|".join(SAMPLES)
    shell:
        "ln -sf $(realpath {input.r1}) {output.r1} && "
        "ln -sf $(realpath {input.r2}) {output.r2}"

##################################################
# Processing

rule cutadapt_primers:
    input:
        r1="fastq/raw/{sid}_R1.fastq.gz",
        r2="fastq/raw/{sid}_R2.fastq.gz"
    output:
        r1="fastq/cutadapt_primers/{sid}_R1.fastq.gz",
        r2="fastq/cutadapt_primers/{sid}_R2.fastq.gz"
    log:
        "fastq/cutadapt_primers/{sid}.log"
    wildcard_constraints:
        sid="|".join(SAMPLES)
    threads:
        config["cutadapt"]["threads"]
    params:
        primer_f=config["primer"]["forward"],
        primer_r_rc=config["primer"]["reverse_revcompl"],
        primer_r=config["primer"]["reverse"],
        primer_f_rc=config["primer"]["forward_revcompl"]
    conda:
        os.path.join(ENV_DIR, "cutadapt.yaml")
    message:
        "Cutadapt: {wildcards.sid}"
    shell:
        # parameters: see wiki for details
        "cutadapt -j {threads} --pair-filter=any {config[cutadapt][params]} "
        "-a '{params.primer_f}...{params.primer_r_rc}' -A '{params.primer_r}...{params.primer_f_rc}' "
        "-o {output.r1} -p {output.r2} {input.r1} {input.r2} &> {log}"

if config["cutadapt"]["strategy"] == "TruSeq":
    rule cutadapt_adapters:
        input:
            r1="fastq/cutadapt_primers/{sid}_R1.fastq.gz",
            r2="fastq/cutadapt_primers/{sid}_R2.fastq.gz"
        output:
            r1="fastq/cutadapt_adapters/{sid}_R1.fastq.gz",
            r2="fastq/cutadapt_adapters/{sid}_R2.fastq.gz"
        log:
            "fastq/cutadapt_adapters/{sid}.log"
        wildcard_constraints:
            sid="|".join(SAMPLES)
        threads:
            config["cutadapt"]["threads"]
        conda:
            os.path.join(ENV_DIR, "cutadapt.yaml")
        message:
            "Cutadapt: {wildcards.sid}"
        shell:
            "cutadapt -j {threads} --pair-filter=any {config[cutadapt][params]} "
            # 3' adapters for R1 (-a) and R2 (-A)
            "-a \"{config[cutadapt][adapters][forward]}\" -A \"{config[cutadapt][adapters][reverse]}\" "
            "-o {output.r1} -p {output.r2} {input.r1} {input.r2} &> {log}"


##################################################
# QC

# FastQC
rule fastqc:
    input:
        "fastq/{ftype}/{sid}_{rid}.fastq.gz"
    output:
        zip="fastqc/{ftype}/{sid}_{rid}_fastqc.zip",
        html="fastqc/{ftype}/{sid}_{rid}_fastqc.html"
    log:
        "fastqc/{ftype}/{sid}_{rid}.log"
    wildcard_constraints:
        ftype="raw|cutadapt_primers|cutadapt_adapters",
        sid="|".join(SAMPLES),
        rid="R1|R2"
    threads:
        config["fastqc"]["threads"]
    conda:
        os.path.join(ENV_DIR, "fastqc.yaml")
    message:
        "FastQC: {input}"
    shell:
        "fastqc -q -f fastq -t {threads} -o $(dirname {output.zip}) {input} &> {log}"

# FastQC summary
rule multiqc_fastqc:
    input:
        expand("fastqc/{{ftype}}/{sid}_{rid}_fastqc.zip", sid=SAMPLES, rid=["R1", "R2"])
    output:
        html="multiqc/fastqc/{ftype}/multiqc_report.html",
        stat="multiqc/fastqc/{ftype}/multiqc_data/multiqc_fastqc.txt",
    log:
        "multiqc/fastqc/{ftype}/multiqc.log"
    wildcard_constraints:
        ftype="raw|cutadapt_primers|cutadapt_adapters",
    threads: 1
    conda:
        os.path.join(ENV_DIR, "multiqc.yaml")
    message:
        "MultiQC (FastQC): {wildcards.ftype}"
    shell:
        "multiqc --interactive -p -f -m fastqc -o $(dirname {output.html}) $(dirname {input[0]}) &> {log}"

# Cutadapt summary
rule multiqc_cutadapt:
    input:
        expand("fastq/{{ftype}}/{sid}_{rid}.fastq.gz", sid=SAMPLES, rid=["R1", "R2"])
    output:
        html="multiqc/fastq/{ftype}/multiqc_report.html"
    log:
        "multiqc/fastq/{ftype}/multiqc.log"
    wildcard_constraints:
        ftype="cutadapt_primers|cutadapt_adapters",
    threads: 1
    conda:
        os.path.join(ENV_DIR, "multiqc.yaml")
    message:
        "MultiQC (Cutadapt)"
    shell:
        "multiqc --interactive -p -f -m cutadapt -o $(dirname {output.html}) $(dirname {input[0]}) &> {log}"
