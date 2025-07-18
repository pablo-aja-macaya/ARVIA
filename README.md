# ARVIA

- [] Herramienta variant calling p. aeruginosa    
    - Nombres
        ARVIA: Antibiotic Resistance Variant Identifier for Pseudomonas aeruginosa
        PAVRA: Pseudomonas aeruginosa Variants and Resistance Analyzer
        PAVCRA: Pseudomonas Aeruginosa Variant Calling Resistance Analysis
        PARVI: P. Aeruginosa Resistance Variant Inspector
    - Funciones:
        - [] Input: paired reads, long reads or assembly
        - [] Output: tabla comparativa a lo ancho (.xlsx y .tsv), tabla comparativa a lo largo (.xlsx y .tsv), informe html de igvvariant, parameter log
        - [] Funciones destacables: 
            - paired reads, long reads or assembly (specify which one was used in each sample)
            - mutations in specific genes related to resistance
            - oprD with closest reference
            - missing features
            - possible polymorphisms
            - truncated genes (only with assembly!)
            - mixed variants
            - comparatives
            - automatic reference download
            - snakemake progress bar (hideable)
            - tests
            - add rgi and mlst?
    - Dependencies:
        - python
        - snakemake
        - snippy
        - bwa
        - samtools
        - makeblastdb y blast (version in bakta)
        - minimap2 for long reads?
        - rgi?
        - mlst?
    - Usage:


## Usage

```sh
# Full pipeline (reads+assemblies)
arvia --input_assemblies folder/*.fasta --input_reads folder/*.fastq.gz --output_folder arvia --threads 60

# Full pipeline using only assemblies (no depth inference in variant calling)
arvia --input_assemblies folder/*.fasta --output_folder arvia --threads 60

# Partial pipeline using only reads (truncation information in assembly from assembly is missing)
arvia --input_reads folder/*.fastq.gz --output_folder arvia --threads 60
```

## Installation

```sh
# Create environment
mamba create -n arvia \
    snakemake==7.18.0 python==3.8.10 pandas==1.5.0 numpy==1.23.1 'biopython>=1.78' rich-argparse==1.6.0 'colorama==0.4.4' 'odfpy==1.4.1' 'setuptools<=70' \
    seqkit==2.1.0 'pigz>=2.4' \
    perl-bioperl snippy==4.6.0 snpEff==4.3.1t bcftools=1.15 openssl==3.5.0 samtools=1.18 blast=2.16.0

git clone https://github.com/Pablo-Aja-Macaya/ARVIA.git
cd ARVIA
python setup.py develop


# gitpython cgecore kma tabulate mlst-cge=2.0.9
```

