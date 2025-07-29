<p align="left">
  <!-- <img src="arvia/data/arvia_header.png" height="70" > -->
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="arvia/data/arvia_header_bb.png">
    <source media="(prefers-color-scheme: light)" srcset="arvia/data/arvia_header_wb.png">
    <img alt="ARVIA Logo" src="xxx" height="70">
  </picture>
</p>

## Summary

ARVIA (**A**ntibiotic **R**esistance **V**ariant **I**dentifier for *Pseudomonas **a**eruginosa*) takes **single-end/paired-end reads (long or short)** and/or an **assembly** per sample to perform exhaustive variant calling of genes related to antibiotic resistance in *Pseudomonas aeruginosa*. See [Usage](#usage) and [Installation](#installation) sections. Its main functions are:
- **Point mutations (SNVs, indels, frameshifts) in PAO1**.
- **Variant calling of closest oprD reference**. 
- Detection of possible **missing features** (e.g. lost genes due to chromosomic rearrangement).
- Detection of possible **truncated genes** due to big chromosomic rearrangements (only with assembly!).
- Detection of **mixed positions** (e.g. 50% of reads indicate C and the other 50% T).
- Detection of possible **polymorphisms** that do not influence antibiotic resistance.
- Creation of **comparative tables** to more easily assess the cause of different phenotypes between samples.


<p align="center">
  <img src="arvia/data/arvia_sp_v0.1.1.png" height="600" >
</p>


## Index
- [Usage](#usage)
- [Installation](#installation)
- [Input](#input)
  - [Input YAML convention](#input-yaml-convention)
  - [ARVIA's naming convention](#file-naming-convention)
- [Full command list](#full-command-list)
- [Citation](#citation) 



## Usage

You can run ARVIA easily with an `input.yaml` file (see [**Input YAML convention**](#input-yaml-convention)) containing the input files:

```sh
# Run ARVIA
arvia run --input_yaml input.yaml --output_folder arvia
```

If your files follow [**ARVIA's naming convention**](#file-naming-convention), you can also give them all with `--reads` and/or `--assemblies` and ARVIA will associate each file to their `sample_id`:

```sh
# Full pipeline (reads+assemblies)
arvia run --assemblies folder/*.fasta --reads folder/*.fastq.gz --output_folder arvia

# Full pipeline using only assemblies (no depth inference in variant calling)
arvia run --assemblies folder/*.fasta --output_folder arvia

# Partial pipeline using only reads (truncation information in assembly from assembly is missing)
arvia run --reads folder/*.fastq.gz --output_folder arvia
```

> [!TIP]
> You can also previsualize what the pipeline is going to do with `--previsualize`:
>```sh
># Run ARVIA
>arvia run --input_yaml input.yaml --output_folder arvia --previsualize
>```

Check out more options, like `--cores`, in the [Full command list](#full-command-list).

## Installation

> [!NOTE]
> This application has been **designed for Linux systems** and tested in **Ubuntu**.

Installation through mamba is highly recommended:

```sh
# Create main environment where ARVIA runs
mamba create -n arvia \
    snakemake==7.18.0 python=3.8 pandas==1.5.0 numpy==1.23.1 'biopython>=1.78' rich-argparse==1.6.0 'colorama==0.4.4' 'odfpy==1.4.1' 'setuptools<=70' toml==0.10.2 xlsxwriter ipykernel \
    seqkit==2.1.0 'pigz>=2.4' ncbi-amrfinderplus mlst unzip \
    perl-bioperl snippy==4.6.0 snpEff==4.3.1t bcftools=1.21 openssl==3.5.0 samtools=1.18 blast=2.16.0
    
conda activate arvia

git clone https://github.com/Pablo-Aja-Macaya/ARVIA.git
cd ARVIA
python -m pip install -e . # "-e" allows for editable mode, else "python -m pip install ."

# Other environments that can't go into the previous one due to incompatibility
mamba create -n arvia_mlst \
    mlst-cge==2.0.9

mamba create -n arvia_rgi \
    rgi=5.2.0



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

## Input

ARVIA takes **single-end/paired-end reads (long or short)** and/or an **assembly** for each sample given. Single-end reads will be considered long reads such as PacBio or Oxford Nanopore Technologies (ONT). It needs at least one of the two types of files, with a maximum of 1 assembly and 2 reads files per sample.

> [!IMPORTANT]
> Selected **pipeline depends on user input**. Every part is available for each input type except the detection of truncated genes caused by big reordenations, which requires an assembly. It is **recommended to provide reads and an assembly** for a more in-depth analysis!


### Input YAML convention

In order to use `--input_yaml` generate a YAML file with the following structure, where keys are unique sample_ids, containing one or two lists named `reads` and `assembly` with their corresponding files (path can be relative from where ARVIA is executed):

```yaml
# -- Input template --
# SAMPLE_10: # <- Will be used as ID
#   reads:
#     - path/to/blablabla_R1.fastq.gz
#     - path/to/blablabla_R2.fastq.gz
#   assembly:
#     - path/to/blablabla.fasta

# -- Your samples --
# Complete example with paired reads and assembly
ARGA00024:
  reads:
    - input/ARGA00024_R2.fastq.gz
    - input/ARGA00024_R1.fastq.gz
  assembly:
    - input/ARGA00024.fasta

# Example with only single-end long reads
# you dont need to specify assembly key if you dont have it
ARGA00461:
  reads:
    - input/ARGA00461.fastq.gz

# Example with only assembly
# you dont need to specify reads key if you dont have it
ARGA00461-a:
  assembly:
    - input/ARGA00461.fasta
```


### File naming convention

You can see the convention expected for `--reads` and `--assemblies` with `--help`:

```sh
-r, --reads path [path ...]         Input reads files. Can be paired-end or single-end and must follow one of these
                                    structures: '{sample_id}.fastq.gz' / '{sample_id}_R[1,2].fastq.gz' /
                                    '{sample_id}_[1,2].fastq.gz' / '{sample_id}_S\d+_L\d+_R[1,2]_\d+.fastq.gz'
-a, --assemblies path [path ...]    Input assembly files. Must follow one of these structures:
                                    '{sample_id}.{fasta,fna,fa,fas}' (default: None)
```



## Full command list 
Full command list available with `arvia --help`:

```sh
Usage:  run [-h] [-i path] [-r path [path ...]] [-a path [path ...]] -o path [-c int] [-p] [--use_conda] [--barcodes str [str ...]]
            [--draw_wf str]

ARVIA: Antibiotic Resistance Variant Identifier for Pseudomonas aeruginosa

Optional Arguments:
  -h, --help                                        show this help message and exit

Input/Output:
  -i, --input_yaml path                             Input files from a YAML. Each key is a sample_id containing two lists of paths with
                                                    keys 'reads' and 'assembly' (default: None)
  -r, --reads path [path ...]                       Input reads files. Can be paired-end or single-end and must follow one of these
                                                    structures: '{sample_id}.fastq.gz' / '{sample_id}_R[1,2].fastq.gz' /
                                                    '{sample_id}_[1,2].fastq.gz' / '{sample_id}_S\d+_L\d+_R[1,2]_\d+.fastq.gz' (default:
                                                    None)
  -a, --assemblies path [path ...]                  Input assembly files. Must follow one of these structures:
                                                    '{sample_id}.{fasta,fna,fa,fas}' (default: None)
  -o, --output_folder path                          Output folder (default: None)

Optional Parameters:
  -c, --cores int                                   Number of cores (default is available cores - 1) (default: 63)
  -p, --previsualize                                Previsualize pipeline to see if everything is as expected (default: False)
  --use_conda                                       If True, use conda environment specified by snakefile (default: False)
  --barcodes str [str ...]                          Space separated list of sample IDs. Only these samples will be processed (default: None)
  --draw_wf str                                     Draw pipeline to this path (PDF (default: None)
```


## Test

In order to test ARVIA's installation, execute the following command:

```sh
arvia test --output_folder test_arvia
```

This command downloads a set of reads and assemblies and tries to run the pipeline.

## Citation

Please cite the database from which PAO1 genome and gene information were retrieved, **[Pseudomonas.com](https://www.pseudomonas.com)**:

Winsor GL, Griffiths EJ, Lo R, Dhillon BK, Shay JA, Brinkman FS (2016). Enhanced annotations and features for comparing thousands of Pseudomonas genomes in the Pseudomonas genome database. Nucleic Acids Res. (2016) doi: 10.1093/nar/gkv1227 (Database issue). Pubmed: 26578582

## Development

### Upload to TestPyPi

```
cd ARVIA/
conda activate twine

# Create build dist
python -m build

# Upload to TestPyPi with twine
twine upload --repository testpypi dist/*

# Now you can pip install
pip install -i https://test.pypi.org/simple/ arvia

```


<!-- 
- [] Herramienta variant calling p. aeruginosa    
    - Funciones:
        - [] Input: paired reads, long reads or assembly
        - [] Output: tabla comparativa a lo ancho (.xlsx y .tsv), tabla comparativa a lo largo (.xlsx y .tsv), informe html de igvvariant, parameter log
        - [] To-do    
            - [] automatic reference download
            - [] in results_per_sample
                - [] format blast table (add header at least)
                - [X] add original muts without filters
            - [] hideable snakemake progress bar?
            - [] tests
            - [] informe html de igvvariant
            - [] rgi
            - [] mlst
            - [] add approximate depth if using reads
            - [] en tabla final si no ha habido ningún cambio en un gen este no aparece, arreglar y meterlo sí o sí aunque esté vacío
            - [] in xlsx output check it looks good on every platform (breaks like \n dont work in windows)
            - [] quitar lo de func en la tabla de parámetros no sé qué es
    - Dependencies:
        - python
        - snakemake
        - snippy
        - bwa
        - samtools
        - makeblastdb y blast (version in bakta)
        - minimap2 for long reads?
        - rgi?
        - mlst-cge?
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

<!-- [home]: https://github.com/Pablo-Aja-Macaya/ARVIA
[installation]: https://github.com/Pablo-Aja-Macaya/ARVIA/#installation



