##################################################
# KRAKEN2

rule KRAKEN2:
    input:
        r1=expand("fastq/cutadapt_adapters/{sid}_R1.fastq.gz", sid=SAMPLES),
        r2=expand("fastq/cutadapt_adapters/{sid}_R2.fastq.gz", sid=SAMPLES),
        database=config['kraken2']['db']
    output:
        rep=os.path.join("kraken2/{sid}.kraken.report.txt"),
        summary=os.path.join("kraken2/{sid}.kraken.summary.out")
    log:
        "kraken2/{sid}.log"
    threads:
        config["kraken2"]["threads"]
    conda:
        os.path.join(ENV_DIR, "kraken2.yaml")
    message:
        "KRAKEN2: {wildcards.sid}"
    shell:
        "(date && kraken2 --threads {threads} --db {input.database} --use-names --confidence 0.5 --paired {input.r1} {input.r2} --gzip-compressed --output {output.summary} --report {output.rep} && date) &> {log}"
