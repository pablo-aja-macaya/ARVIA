import pandas as pd
import glob
import json
import re
import traceback
from pathlib import Path
from snakemake.logging import logger
import logging
import warnings
import datetime
import arvia

from arvia.utils.aeruginosa_snippy import filter_snippy_result
from arvia.utils.console_log import CONSOLE_STDOUT, CONSOLE_STDERR, log_error_and_raise
from arvia.utils.annotation_extraction import get_proteins_from_gbk

from arvia.utils.local_paths import OPRD_NUCL, OPRD_CONFIG #, KPNEUMONIAE_SELECTED_GENES
from arvia.utils.local_paths import PAERUGINOSA_GENOME_GBK
from arvia.utils.local_paths import CONDA_ENVS
# from bactasys.utils.config.local_paths import MLST_DB, MLST_CONFIG
# from bactasys.utils.config.local_paths import CARD_JSON

warnings.simplefilter(action='ignore', category=FutureWarning) # remove warning from pandas
ARVIA_DIR = arvia.__file__.replace("/__init__.py", "")  # get install directory of bactasys
DATETIME_OF_CALL = datetime.datetime.now()



if config:
    snakemake_console_log = config.get("snakemake_console_log")

    # # ---- Get parameters ----
    # # Input
    # SHORT_READS_INPUT = config["sr_folder"]
    # LONG_READS_INPUT = config["lr_folder"]
    # # READ_TYPE = ? # TODO: this

    # # Output
    # BASE_FOLDER = config["output_folder"]

    # # Filter barcodes
    # USER_BARCODES = config.get("barcodes")

# ---- Send snakemake log to custom file to parse with custom log ----
if snakemake_console_log is not None:
    handler = logging.FileHandler(snakemake_console_log, mode='a')
    formatter = logging.Formatter('%(message)s')
    handler.setFormatter(formatter)
    logger.set_stream_handler(handler)


# ---- Input set-up ----
temp_input_assembly_folder = "/home/usuario/Proyectos/Results/ARGA/ARGA_ALL_v2/assembly/05_final_assembly"
temp_input_shortreads_folder = "/home/usuario/Proyectos/Results/ARGA/ARGA_ALL_v2/short_reads_qc/06_clean_reads"
temp_input_longreads_folder = "/home/usuario/Proyectos/Results/ARGA/ARGA_ALL_v2/00_long_clean_reads"
INPUT_FILES = {
    # Single reads with or without assembly
    "ARGA00461": { # full pipeline
        "reads": [f"{temp_input_longreads_folder}/ARGA00461/ARGA00461.fastq.gz"],
        "assembly": [f"{temp_input_assembly_folder}/ARGA00461/ARGA00461_assembly.fasta"]
    },
    "ARGA00190": { # truncation cant be done
        "reads": [f"{temp_input_longreads_folder}/ARGA00190/ARGA00190.fastq.gz"],
        "assembly": []
    },
    # Paired reads with or without assembly
    "ARGA00024": { # full pipeline
        "reads": [f"{temp_input_shortreads_folder}/ARGA00024/ARGA00024_R1.fastq.gz", f"{temp_input_shortreads_folder}/ARGA00024/ARGA00024_R2.fastq.gz"],
        "assembly": [f"{temp_input_assembly_folder}/ARGA00024/ARGA00024_assembly.fasta"]
    },
    "ARGA00025": { # truncation cant be done
        "reads": [f"{temp_input_shortreads_folder}/ARGA00025/ARGA00025_R1.fastq.gz", f"{temp_input_shortreads_folder}/ARGA00025/ARGA00025_R2.fastq.gz"],
        "assembly": [f"{temp_input_assembly_folder}/ARGA00025/ARGA00025_assembly.fasta"]
    },
    # Assembly without reads
    "ARGA00031": { # everything with assembly
        "reads": [],
        "assembly": [f"{temp_input_assembly_folder}/ARGA00031/ARGA00031_assembly.fasta"]
    },
}
  
for k,v in INPUT_FILES.items():
    reads = v["reads"]
    assembly = v["assembly"]

    assert type(reads)==list, "Reads are not in list format"
    assert type(assembly)==list, "Assembly is not in list format"
    check_reads_exists = [i for i in reads if not Path(i).exists()]
    check_assembly_exists = [i for i in assembly if not Path(i).exists()]

    # If reads are supplied
    if (len(reads) in [1,2]):
        # assert len(check_reads_exists)==0, f"At least a file path does not exist: {check_reads_exists}"
        v["reads_type"] = "single_end" if len(reads)==1 else "paired_end"

        # If assembly is supplied
        if assembly:
            # assert check_assembly_exists, f"At least a file path does not exist: {assembly}"
            v["pipeline"] = "full_pipeline"
        else:
            v["pipeline"] = "only_reads"
    
    elif assembly:
        # assert check_assembly_exists, f"At least a file path does not exist: {assembly}"
        v["pipeline"] = "only_assembly"
        v["reads_type"] = None
    
    else:
        raise Exception(f"Unexpected input -> {k}: {v}")

# CONSOLE_STDOUT.print(INPUT_FILES)



# ---- Output folders ----
PIPELINE_OUTPUT = config["output_folder"]

# Input
GET_ASSEMBLIES_OUTPUT = f"{PIPELINE_OUTPUT}/00_input_assemblies"
GET_READS_OUTPUT = f"{PIPELINE_OUTPUT}/00_input_reads"

# General snippy
PAERUGINOSA_MUTS_OUTPUT = f"{PIPELINE_OUTPUT}/variant_calling/run"
PAERUGINOSA_MUTS_PROCESS_OUTPUT = f"{PIPELINE_OUTPUT}/variant_calling/process"

# OprD
PAERUGINOSA_OPRD_OUTPUT = f"{PIPELINE_OUTPUT}/oprd"
ALIGN_OPRD = f"{PAERUGINOSA_OPRD_OUTPUT}/01_align"
DECIDE_BEST_OPRD = f"{PAERUGINOSA_OPRD_OUTPUT}/02_decide"
SNIPPY_OPRD = f"{PAERUGINOSA_OPRD_OUTPUT}/03_snippy"

# Truncated genes
PAERUGINOSA_TRUNC_OUTPUT = f"{PIPELINE_OUTPUT}/truncations"
EXTRACT_PAERUGINOSA_REF_GENES_OUTPUT = f"{PAERUGINOSA_TRUNC_OUTPUT}/01_extract_ref_genes"
MAKEBLASTDB_FROM_ASSEMBLY_OUTPUT = f"{PAERUGINOSA_TRUNC_OUTPUT}/02_blastdb"
BLAST_PAERUGINOSA_GENES_TO_ASSEMBLY_OUTPUT = f"{PAERUGINOSA_TRUNC_OUTPUT}/03_blast"

RESULTS_PER_SAMPLE = f"{PIPELINE_OUTPUT}/reunite"

# ---- Params ----


##########################################
# --------------- Rules ---------------- #
##########################################


# ---- Common ----
def get_if_use_assembly_or_reads(wc):
    bc_reads_type = INPUT_FILES[wc.barcode]["reads_type"]
    return bc_reads_type if bc_reads_type else "assembly"

def get_input_assembly(wc):
    return INPUT_FILES[wc.barcode]["assembly"]

def get_input_reads(wc):
    return INPUT_FILES[wc.barcode]["reads"]

rule snippy:
    input:
        assembly="",
        reads=[],
        ref_gbk=""
    output:
        folder=directory(Path("folder")),
        res=Path(".tab"),
        gene_coverage=Path(".tsv"),
    params:
        selected_input=None, # "paired_end" | "single_end" | "assembly",
        min_depth=config.get("snippy", {}).get("min_depth", 5),
        maxsoft=config.get("snippy", {}).get("maxsoft", 1000),
        arvia_dir=ARVIA_DIR,
    threads: 6
    conda:
        CONDA_ENVS["arvia"]
    log:
        Path("", "arvia.log")
    shell:
        """
        (
        rm -r {output.folder}

        if [[ {params.selected_input} == "paired_end" ]]; then
            reads=({input.reads})
            snippy --R1 ${{reads[0]}} --R2 ${{reads[1]}} --ref {input.ref_gbk} --mincov {params.min_depth} --maxsoft {params.maxsoft} --outdir {output.folder} --cpus {threads} --quiet

        elif [[ {params.selected_input} == "single_end" ]]; then
            snippy --se {input.reads} --ref {input.ref_gbk} --mincov {params.min_depth} --maxsoft {params.maxsoft} --outdir {output.folder} --cpus {threads} --quiet

        elif [[ {params.selected_input} == "assembly" ]]; then
            snippy --ctgs {input.assembly} --ref {input.ref_gbk} --mincov {params.min_depth} --maxsoft {params.maxsoft} --outdir {output.folder} --cpus {threads} --quiet
        else
            exit 1
        fi

        # Run fixed vcf-to-tab script
        {params.arvia_dir}/scripts/snippy-vcf_to_tab --ref {output.folder}/reference/ref.fa --gff {output.folder}/reference/ref.gff --vcf {output.folder}/snps.vcf > {output.folder}/snps.tab

        # ---- Get mixed variants (no filter) ----
        # Difference: quitar el FMT/GT="1/1" del bcftools view
        cd {output.folder}
        bcftools view --include 'QUAL>=100 && FMT/DP>=5 && (FMT/AO)/(FMT/DP)>=0.2' snps.raw.vcf  | vt normalize -r reference/ref.fa - | bcftools annotate --remove '^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL' > snps.nofilt.vcf 2>> snps.log
        snpEff ann -noLog -noStats -no-downstream -no-upstream -no-utr -c reference/snpeff.config -dataDir . ref snps.nofilt.vcf > snps.nofilt.final.vcf 2>> snps.log
        {params.arvia_dir}/scripts/snippy-vcf_to_tab --gff reference/ref.gff --ref reference/ref.fa --vcf snps.nofilt.final.vcf > snps.nofilt.tab 2>> snps.log

        # ---- Get gene coverage ----
        cd {output.folder}
        bash {params.arvia_dir}/scripts/snippy_add_missing_features.sh

        # ---- Clean up ----
        rm {output.folder}/reference -r
        rm {output.folder}/*.bam {output.folder}/*.bai {output.folder}/*.html {output.folder}/*.gff {output.folder}/snps.subs.vcf {output.folder}/ref.fa.fai {output.folder}/*.fa
        ) &> {log}
        """




# ---- Pseudomonas aeruginosa Snippy PAO1 ----
use rule snippy as paeruginosa_mutations with:
    input:
        assembly=lambda wc: get_input_assembly(wc),
        reads=lambda wc: get_input_reads(wc),
        ref_gbk=PAERUGINOSA_GENOME_GBK,
    output:
        folder=directory(Path(PAERUGINOSA_MUTS_OUTPUT, "{barcode}")),
        res=Path(PAERUGINOSA_MUTS_OUTPUT, "{barcode}", "snps.tab"),
        gene_coverage=Path(PAERUGINOSA_MUTS_OUTPUT, "{barcode}", "gene_coverage.tsv"),
    params:
        selected_input=lambda wc: get_if_use_assembly_or_reads(wc), # "paired_end" | "single_end" | "assembly",
        min_depth=config.get("snippy", {}).get("min_depth", 5),
        maxsoft=config.get("snippy", {}).get("maxsoft", 1000),
        arvia_dir=ARVIA_DIR,
    threads: 6
    log:
        Path(PAERUGINOSA_MUTS_OUTPUT, "{barcode}", "arvia.log")

rule process_paeruginosa_mutations:
    input:
        res=rules.paeruginosa_mutations.output.res,
    output:
        folder=directory(Path(PAERUGINOSA_MUTS_PROCESS_OUTPUT, "{barcode}")),
        res=Path(PAERUGINOSA_MUTS_PROCESS_OUTPUT, "{barcode}", "{barcode}_filtered_snps.tab"),
    threads: 1
    run:
        df = filter_snippy_result(input.res, output.res)


# ---- Pseudomonas aeruginosa blast reference genes ----
rule extract_paeruginosa_ref_genes:
    input:
        ref_gbk=PAERUGINOSA_GENOME_GBK,
    output:
        folder=directory(Path(EXTRACT_PAERUGINOSA_REF_GENES_OUTPUT, "genes")),
        genes_ffn=Path(EXTRACT_PAERUGINOSA_REF_GENES_OUTPUT, "genes","genes.ffn"),
        genes_faa=Path(EXTRACT_PAERUGINOSA_REF_GENES_OUTPUT, "genes","genes.faa"),
    run:
        _ = get_proteins_from_gbk(input.ref_gbk, output.genes_ffn, output_aa=False)
        _ = get_proteins_from_gbk(input.ref_gbk, output.genes_faa, output_aa=True)


rule makeblastdb_from_assembly:
    input:
        ref=lambda wc: get_input_assembly(wc),
    output:
        folder=directory(Path(MAKEBLASTDB_FROM_ASSEMBLY_OUTPUT, "{barcode}")),
    params:
        blast_db="db",
        dbtype="nucl",  # nucl | prot
    conda:
        CONDA_ENVS["arvia"] # uses this but it just needs blast
    threads: 20
    log:
        Path(MAKEBLASTDB_FROM_ASSEMBLY_OUTPUT, "{barcode}", "arvia.log")
    shell:
        """
        mkdir -p {output.folder}
        (
        cat {input.ref} | makeblastdb -dbtype {params.dbtype} -parse_seqids -title {output.folder} -out {output.folder}/{params.blast_db}
        ) &> {log}
        """


rule blast_paeruginosa_genes_to_assembly:
    input:
        query=rules.extract_paeruginosa_ref_genes.output.genes_ffn,
        db=rules.makeblastdb_from_assembly.output.folder,
    output:
        folder=directory(Path(BLAST_PAERUGINOSA_GENES_TO_ASSEMBLY_OUTPUT,"{barcode}")),
        res=Path(BLAST_PAERUGINOSA_GENES_TO_ASSEMBLY_OUTPUT,"{barcode}","{barcode}.tsv"),
    threads: 20
    conda:
        CONDA_ENVS["arvia"] # uses this but it just needs blast
    params:
        outfmt="6 qseqid sseqid stitle pident qcovs qcovhsp qcovus length mismatch gapopen qlen slen qstart qend sstart send evalue bitscore qseq sseq", # warning: different to the default blst_outfmt used in other rules
        blast_type="blastn",
        max_target_seqs=5,
    log:
        Path(BLAST_PAERUGINOSA_GENES_TO_ASSEMBLY_OUTPUT, "{barcode}", "arvia.log")
    shell:
        """
        mkdir -p {output.folder}
        (
        {params.blast_type} -query {input.query} -out {output.res} -max_target_seqs {params.max_target_seqs} -db {input.db}/db -subject_besthit -num_threads {threads} -outfmt "{params.outfmt}"
        ) &> {log}
        """

# ---- Pseudomonas aeruginosa Snippy Closest oprD ----
rule align_oprd:
    input:
        ref = OPRD_NUCL,
        assembly=lambda wc: get_input_assembly(wc),
        reads=lambda wc: get_input_reads(wc),
    output:
        folder = directory(Path(ALIGN_OPRD, "{barcode}")),
        bam = temp(Path(ALIGN_OPRD, "{barcode}", "aln.bam")),
        coverage = Path(ALIGN_OPRD, "{barcode}", "coverage.tsv"),
    params:
        selected_input=lambda wc: get_if_use_assembly_or_reads(wc), # "paired_end" | "single_end" | "assembly",
    threads: 12
    conda:
        CONDA_ENVS["arvia"]
    log:
        Path(ALIGN_OPRD, "{barcode}", "arvia.log")
    shell:
        """
        mkdir -p {output.folder}

        (
        # Align
        cp {input.ref} {output.folder}/ref.fasta


        if [[ {params.selected_input} == "paired_end" ]]; then
            bwa index {output.folder}/ref.fasta;     
            bwa mem -t {threads} {output.folder}/ref.fasta {input.reads} | samtools view --threads {threads} -b - | samtools sort --threads {threads} - > {output.bam} 

        elif [[ {params.selected_input} == "single_end" ]]; then
            minimap2 -a -x map-hifi -t {threads} {output.folder}/ref.fasta {input.reads} | samtools view --threads {threads} -b - | samtools sort --threads {threads} - > {output.bam} 

        elif [[ {params.selected_input} == "assembly" ]]; then
            minimap2 -a -x asm10 -t {threads} {output.folder}/ref.fasta {input.assembly} | samtools view --threads {threads} -b - | samtools sort --threads {threads} - > {output.bam} 

        else
            exit 1
        fi

        samtools index {output.bam} 

        # Decide best reference
        samtools coverage {output.bam} > {output.coverage}

        # Clean 
        rm {output.folder}/ref.fast* {output.folder}/*.bai
        ) &> {log}
        """

rule decide_best_oprd_ref:
    input:
        coverage = rules.align_oprd.output.coverage,
    output:
        folder = directory(Path(DECIDE_BEST_OPRD, "{barcode}")),
        selected_ref = temp(Path(DECIDE_BEST_OPRD, "{barcode}", "ref.gbk")),
        selected_ref_txt = Path(DECIDE_BEST_OPRD, "{barcode}", "selected_ref.txt"),
    params:
        oprd_config = OPRD_CONFIG
    run:
        df = pd.read_csv(input.coverage, sep="\t")
        df = df.sort_values(["coverage","meandepth"], ascending=False)

        locus_tag = df.iloc[0]["#rname"]
        coverage = df.iloc[0]["coverage"]
        meandepth = df.iloc[0]["meandepth"]

        strain = params.oprd_config["oprD"][locus_tag]["strain"]
        gbk = params.oprd_config["oprD"][locus_tag]["gbk"]

        # CONSOLE_STDOUT.log(f"Best match is {locus_tag} ({strain}) at {coverage}% and {int(meandepth)}x. GBK: {gbk}")
        shell("cp {gbk} {output.selected_ref}")

        with open(output.selected_ref_txt, "wt") as handle:
            handle.write(f"{wildcards.barcode}\t{locus_tag}\t{strain}\t{coverage}\t{int(meandepth)}\n")

use rule snippy as paeruginosa_oprd with:
    input:
        assembly=lambda wc: get_input_assembly(wc),
        reads=lambda wc: get_input_reads(wc),
        ref_gbk=rules.decide_best_oprd_ref.output.selected_ref,
    output:
        folder = directory(Path(SNIPPY_OPRD, "{barcode}")),
        res = Path(SNIPPY_OPRD, "{barcode}", "snps.tab"),
        gene_coverage=Path(SNIPPY_OPRD, "{barcode}", "gene_coverage.tsv"),
    params:
        selected_input=lambda wc: get_if_use_assembly_or_reads(wc), # "paired_end" | "single_end" | "assembly",
        min_depth=config.get("snippy", {}).get("min_depth", 5),
        maxsoft=config.get("snippy", {}).get("maxsoft", 1000),
        arvia_dir=ARVIA_DIR,
    threads: 12
    log:
        Path(SNIPPY_OPRD, "{barcode}", "arvia.log")



# ---------------
def decide_steps(wc):
    # use_assembly_or_reads = get_if_use_assembly_or_reads(wc)
    # bc_reads_type = INPUT_FILES[wc.barcode]["reads_type"]
    pipeline = INPUT_FILES[wc.barcode]["pipeline"] # full_pipeline | only_reads | only_assembly
    if pipeline in ["full_pipeline", "only_assembly"]:
        steps = {
            "process_paeruginosa_mutations": rules.process_paeruginosa_mutations.output.res,
            "paeruginosa_gene_coverage": rules.paeruginosa_mutations.output.gene_coverage,
            "paeruginosa_oprd": rules.paeruginosa_oprd.output.res,
            "paeruginosa_oprd_refs": rules.decide_best_oprd_ref.output.selected_ref_txt,
            "paeruginosa_assembly_blast": rules.blast_paeruginosa_genes_to_assembly.output.res,
        }
    elif pipeline == "only_reads":
        steps = {
            "process_paeruginosa_mutations": rules.process_paeruginosa_mutations.output.res,
            "paeruginosa_gene_coverage": rules.paeruginosa_mutations.output.gene_coverage,
            "paeruginosa_oprd": rules.paeruginosa_oprd.output.res,
            "paeruginosa_oprd_refs": rules.decide_best_oprd_ref.output.selected_ref_txt,
            # "paeruginosa_assembly_blast": NULL,
        }
    else:
        raise Exception

    # CONSOLE_STDOUT.print(steps)
    return steps

rule get_results_per_sample:
    input:
        unpack(lambda wc: decide_steps(wc))
    output:
        folder = directory(Path(RESULTS_PER_SAMPLE, "{barcode}")),
    shell:
        """
        mkdir -p {output.folder}
        cp {input} {output.folder}
        """


rule all:
    input:
        expand(rules.get_results_per_sample.output.folder, zip, barcode=list(INPUT_FILES.keys()))
        # rules.merge_species_specific_results.output.folder,
    default_target: True


# ---- Merge ----
# Full or assembly
# Only reads

# rule merge_species_specific_results:
#     input:
#         json_files = expand(rules.get_species_specific_results.output.input_files, **wildcards_dict),
#     output:
#         folder = directory(Path(MERGE_SPECIES_SPECIFIC_RESULTS)),
#     run:
#         shell("mkdir -p {output.folder}")

#         # Read JSONs into list of dicts
#         l = []
#         for i in input.json_files:
#             with open(i) as json_handle:
#                 d = json.load(json_handle)
#                 l.append(d)
        
#         # Combine values if their keys match
#         # Example:
#         #   Input: [{"a": 1, "b": 2},{"a": 5, "b": 7, "c": 0,}] 
#         #   Output: {'a': [1, 5], 'b': [2, 7], 'c': [0]}
#         combined_jsons = {}
#         for d in l:
#             for k,v in d.items():
#                 if combined_jsons.get(k) is None:
#                     combined_jsons[k] = [v]
#                 else:
#                     combined_jsons[k] += [v]

#         # Paeruginosa merge (from annotation.smk's SPECIES_SPECIFIC_ANALYSES)
#         paeruginosa_expected_keys = [
#             "process_paeruginosa_mutations", "paeruginosa_gene_coverage",
#             "paeruginosa_oprd", "paeruginosa_oprd_refs", "paeruginosa_assembly_blast",
#         ]
#         paeruginosa_results_exist = len([i for i in paeruginosa_expected_keys if i not in combined_jsons])==0
#         if paeruginosa_results_exist:
#             from bactasys.utils.genomics.combine_snippy_results import get_default_snippy_combination, paeruginosa_combine_all_mutational_res

#             default_df, bcs = get_default_snippy_combination(
#                 combined_jsons["process_paeruginosa_mutations"], 
#                 output.folder, 
#                 selected_bcs=[], 
#                 paeruginosa=True
#             )
#             paeruginosa_combine_all_mutational_res(
#                 default_df=default_df,
#                 oprd_fs=combined_jsons["paeruginosa_oprd"],
#                 oprd_refs_fs=combined_jsons["paeruginosa_oprd_refs"],
#                 gene_coverage_fs=combined_jsons["paeruginosa_gene_coverage"],
#                 truncation_fs=combined_jsons["paeruginosa_assembly_blast"],
#                 bcs=bcs,
#                 output_folder=output.folder,
#                 filter_poly=False,
#             )

