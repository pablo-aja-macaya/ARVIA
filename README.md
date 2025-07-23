<p align="left">
  <img src="arvia/data/arvia_icon2.png" height="80" >
</p>

## Summary

ARVIA (**A**ntibiotic **R**esistance **V**ariant **I**dentifier for *Pseudomonas **a**eruginosa*) takes **single-end/paired-end reads** and/or **assemblies** to perform exhaustive variant calling of genes related to antibiotic resistance in *Pseudomonas aeruginosa*. Its main functions are:
- **Point mutations (SNVs, indels, frameshifts) in PAO1**.
- **Variant calling of closest oprD reference**. 
- Detection of possible **missing features** (e.g. lost genes due to chromosomic rearrangement).
- Detection of possible **truncated genes** due to big chromosomic rearrangements (only with assembly!).
- Detection of **mixed positions** (e.g. 50% of reads indicate C and the other 50% T).
- Detection of possible **polymorphisms** that do not influence antibiotic resistance.
- Creation of comparative tables.

## Usage

```sh
# Full pipeline (reads+assemblies)
arvia run --assemblies folder/*.fasta --reads folder/*.fastq.gz --output_folder arvia --threads 60

# Full pipeline using only assemblies (no depth inference in variant calling)
arvia run --assemblies folder/*.fasta --output_folder arvia --threads 60

# Partial pipeline using only reads (truncation information in assembly from assembly is missing)
arvia run --reads folder/*.fastq.gz --output_folder arvia --threads 60
```

## Installation

```sh
# Create environment
mamba create -n arvia \
    snakemake==7.18.0 python=3.8 pandas==1.5.0 numpy==1.23.1 'biopython>=1.78' rich-argparse==1.6.0 'colorama==0.4.4' 'odfpy==1.4.1' 'setuptools<=70' toml==0.10.2 xlsxwriter \
    seqkit==2.1.0 'pigz>=2.4' \
    perl-bioperl snippy==4.6.0 snpEff==4.3.1t bcftools=1.15 openssl==3.5.0 samtools=1.18 blast=2.16.0
    
conda activate arvia

git clone https://github.com/Pablo-Aja-Macaya/ARVIA.git
cd ARVIA
python -m pip install -e . # "-e" allows for editable mode, else "python -m pip install ."

```

<!-- 
# Testing package updates (this one works)
mamba create -n arvia_test_env \
    'snakemake==7.18.0' 'python>=3.8.10' 'pandas>=1.5.0' 'numpy>=1.23.1' 'biopython>=1.78' 'rich-argparse>=1.6.0' 'colorama>=0.4.4' 'odfpy>=1.4.1' 'setuptools<81' xlsxwriter \
    seqkit==2.1.0 'pigz>=2.4' \
    perl-bioperl snippy==4.6.0 snpEff==4.3.1t bcftools=1.21 openssl==3.5.0 samtools=1.18 blast=2.16.0 

mamba create -n arvia_test_env \
    'snakemake==9.8.1' 'python>=3.8.10' 'pandas>=1.5.0' 'numpy>=1.23.1' 'biopython>=1.78' 'rich-argparse>=1.6.0' 'colorama>=0.4.4' 'odfpy>=1.4.1' 'setuptools<81' xlsxwriter \
    seqkit==2.1.0 'pigz>=2.4' \
    perl-bioperl snippy==4.6.0 snpEff==4.3.1t bcftools=1.21 openssl==3.5.0 samtools=1.18 blast=2.16.0 

conda activate arvia_test_env
python setup.py develop
-->

## Citation

Please cite the database from which PAO1 genome and gene information were retrieved, **[Pseudomonas.com](https://www.pseudomonas.com)**:

Winsor GL, Griffiths EJ, Lo R, Dhillon BK, Shay JA, Brinkman FS (2016). Enhanced annotations and features for comparing thousands of Pseudomonas genomes in the Pseudomonas genome database. Nucleic Acids Res. (2016) doi: 10.1093/nar/gkv1227 (Database issue). Pubmed: 26578582


<!-- 
- [] Herramienta variant calling p. aeruginosa    
    - Funciones:
        - [] Input: paired reads, long reads or assembly
        - [] Output: tabla comparativa a lo ancho (.xlsx y .tsv), tabla comparativa a lo largo (.xlsx y .tsv), informe html de igvvariant, parameter log
        - [] To-do    
            - [] automatic reference download
            - [] in results_per_sample
                - [] format blast table (add header at least)
                - [] add original muts without filters
            - [] hideable snakemake progress bar
            - [] tests
            - [] add approximate depth if using reads
            - [] revisar parte de blast porque hay genes que no aparecen en tabla final
            - [] in xlsx output check it looks good on every platform (breaks like \n dont work in windows)
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
    - Paper: https://academic.oup.com/bioinformatics/pages/instructions_for_authors
        - [] paper
        - [] cover letter
        - [] Title page
        - [] .tif files (1200 d.p.i. for line drawings and 350 d.p.i. for colour and half-tone artwork). For online submission, please also prepare a second version of your figures at low-resolution for use in the review process; these versions of the figures can be saved in .jpg, .gif, .tif or .eps 
    - Nombres
        ARVIA: Antibiotic Resistance Variant Identifier for Pseudomonas aeruginosa
        PAVRA: Pseudomonas aeruginosa Variants and Resistance Analyzer
        PAVCRA: Pseudomonas Aeruginosa Variant Calling Resistance Analysis
        PARVI: P. Aeruginosa Resistance Variant Inspector
-->
