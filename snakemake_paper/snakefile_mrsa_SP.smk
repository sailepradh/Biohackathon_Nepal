from snakemake.utils import min_version
min_version("5.3.0")

configfile: "config.yml"

rule all:
    """
    Collect the main outputs of the workflow.
    """
    input:
        "results/tables/counts.tsv",
        "results/multiqc.html",
        "results/rulegraph.png"

rule get_SRA_by_accession:
    """
    Retrieve a single-read FASTQ file from SRA (Sequence Read Archive) by run accession number.
    """
    output:
        "data/raw_internal/{sra_id}.fastq.gz"
    params:
        max_reads = config ["max_reads"]
    shell:
        """
        fastq-dump {wildcards.sra_id} -X {params.max_reads} --readids \
            --dumpbase --skip-technical --gzip -Z > {output}
        """

rule fastqc:
    """
    Run FastQC on a FASTQ file.
    """
    input:
        "data/raw_internal/{id}.fastq.gz"
    output:
        "results/{id}_fastqc.html",
        "intermediate/{id}_fastqc.zip"
    shell:
        """
        # Run fastQC and save the output to the current directory
        fastqc {input} -q -o .

        # Move the files which are used in the workflow
        mv {wildcards.id}_fastqc.html {output[0]}
        mv {wildcards.id}_fastqc.zip {output[1]}
        """

rule multiqc:
    """
    Aggregate all FastQC reports into a MultiQC report.
    """
    input:
        expand("intermediate/{sra_id}_fastqc.zip", sra_id =config ["sample_ids"])

    output:
        html = "results/multiqc.html",
        stats = "intermediate/multiqc_general_stats.txt"
    log:
        "results/logs/multiqc/multiqc.log"
#    shadow: "minimal"
    shell:
        """
        # Run multiQC and keep the html report
        multiqc -n multiqc.html {input} 2> {log}
        mv multiqc.html {output.html}
        mv multiqc_data/multiqc_general_stats.txt {output.stats}

        # Remove the other directory that multiQC creates
        rm -rf multiqc_data
        """

rule get_genome_fasta:
    """
    Retrieve the sequence in fasta format for a genome.
    """
    output:
        "data/raw_external/{genome_id}.fa.gz"
    log: 
        "results/log/get_genome_fasta/{genome_id}.log"
    params:
        fasta_path = lambda wildcards: config["genomes"][wildcards.genome_id]["fasta"]
    shell:
        """
        wget {params.fasta_path} -O {output} -o {log}
        """

rule get_genome_gff3:
    """
    Retrieve annotation in gff3 format for a genome.
    """
    output:
        "data/raw_external/{genome_id}.gff3.gz"
    log:
        "results/log/get_genome_gff/{genome_id}.log"
    params:
        gff_path = lambda wildcards: config["genomes"][wildcards.genome_id]["gff3"]
    shell:
        """
        wget {params.gff_path} -O {output} -o {log}
        """

rule index_genome:
    """
    Index a genome using Bowtie 2.
    """
    input:
        "data/raw_external/{genome_id}.fa.gz"
    output:
        expand("intermediate/{{genome_id}}.{my_substr}.bt2", my_substr =  ["1", "2", "3", "4", "rev.1", "rev.2"])
    log:
        "results/logs/index_genome/{genome_id}.log"
 
    shell:
        """
        # Bowtie2 cannot use .gz, so unzip to a temporary file first
        gunzip -c {input} > tempfile
        bowtie2-build tempfile intermediate/{wildcards.genome_id} 2> {log}

        # Remove the temporary file
        rm tempfile
        """

rule align_to_genome:
    """
    Align a fastq file to a genome index using Bowtie 2.
    """
    input:
        "data/raw_internal/{sra_id}.fastq.gz",
        expand("intermediate/{genome_id}.{my_substr}.bt2",
                genome_id = config["genome_id"],
                my_substr = ["1", "2", "3", "4", "rev.1", "rev.2"])
    output:
        temp("intermediate/{sra_id,\w+}.bam")
    log:
        "results/logs/align_to_genome/{sra_id}.log"
    shell:
        """
        genomeid="$(echo {input[1]} | cut -d "." -f1)"
        bowtie2 -x $genomeid -U {input[0]} > {output} 2> {log}
        """

rule sort_bam:
    """
    Sort a bam file.
    """
    input:
        "intermediate/{sra_id}.bam"
    output:
        "intermediate/{sra_id}.sorted.bam"
    shell:
        """
        samtools sort {input} > {output}
        """

rule generate_count_table:
    """
    Generate a count table using htseq-count.
    """
    input:
        bams = expand("intermediate/{sra_id}.sorted.bam", sra_id = config ["sample_ids"]),
        annotation=expand("data/raw_external/{genome_id}.gff3.gz", genome_id = config["genome_id"])
    output:
        "results/tables/counts.tsv"
   # shadow: "minimal"
    shell:
        """
        featureCounts -t exon -g Name -a {input.annotation} -o {output} {input.bams}
        """

rule generate_rulegraph:
    """
    Generate a rulegraph for the workflow.
    """
    output:
        "results/rulegraph.png"
    shell:
        """
        snakemake --snakefile snakefile_mrsa.smk --config max_reads=0 --rulegraph | dot -Tpng > {output}
        """
