import yaml
import re
from pathlib import Path
import pandas as pd


# Possible input file patterns so we can associate the sample name 
INPUT_FILE_PATTERNS = {
    "reads": [
        r"(.+)_S\d+_L\d+_R[12]_\d+.fastq.gz$",
        r"(.+)_R[12].fastq.gz$",
        r"(.+)_[12].fastq.gz$",
        r"(.+).fastq.gz$",
    ],
    "assembly": [
        r"(.+).fasta$",
        r"(.+).fna$",
        r"(.+).fas$",
        r"(.+).fa$",
    ],

}

def associate_user_input_files(config: dict) -> dict:
    """    
    Transform user input (dict) into expected dictionary associating sample ids in name file.
    Sample IDs are extracted from files through a set of possible patterns
    Expected structure below
    """
    d = {
        # "ARGA00461": { 
        #     "reads": ["ARGA00461.fastq.gz"],
        #     "assembly": ["ARGA00461.fasta"]
        # },
    }
    for f in sorted(config["reads"]) + sorted(config["assemblies"]):
        pattern_detected = False
        for pat_type, pats in INPUT_FILE_PATTERNS.items():
            for pat in pats:
                res = re.findall(pat, str(Path(f).name))
                if res:
                    assert len(res)==1 and type(res)==list, "Unexpected"
                    bc = res[0]
                    
                    if not d.get(bc):
                        d[bc] = {
                            "reads": [],
                            "assembly": [],
                        }
                    if pat_type=="reads":
                        d[bc]["reads"] += [f]
                    elif pat_type=="assembly":
                        d[bc]["assembly"] += [f]
                    pattern_detected = True
                    break

        assert pattern_detected, f"Could not find expected file structure in {f}"
    return d

def check_input_file_dict_and_decide_pipeline(d):
    for k,v in d.items():
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
    return d

def expand_input_file_dict_into_multiple_cols(row):
    """
    Takes column 1 and expands it into columns reads_1, reads_2 and assembly
    depending on the number of files and their presence
    """
    read_file_count = len(row[1]["reads"])
    assembly_file_count = len(row[1]["assembly"])

    if read_file_count==2:
        row["reads_1"] = row[1]["reads"][0]
        row["reads_2"] = row[1]["reads"][1]
    elif read_file_count==1:
        row["reads_1"] = row[1]["reads"][0]
    elif read_file_count>2:
        raise Exception(f"Unexpected number of reads: {row}")

    if assembly_file_count == 1:
        row["assembly"] = row[1]["assembly"][0]
    return row

def input_files_dict_to_df(d: dict) -> pd.DataFrame:
    file_manifest_df = pd.DataFrame(d.items())
    file_manifest_df.columns = ["sample", 1]
    file_manifest_df = file_manifest_df.apply(expand_input_file_dict_into_multiple_cols, axis=1)
    file_manifest_df = file_manifest_df.drop(columns=[1])
    file_manifest_df = file_manifest_df.fillna("-")
    return file_manifest_df

# # Read config (user input)
# yaml_config = "/home/usuario/Proyectos/Results/tests/arvia/config.yaml"
# with open(yaml_config, 'r') as f:
#     config = yaml.load(f, Loader=yaml.SafeLoader)

# input_files_dict_to_df(associate_user_input_files(config))


# INPUT_FILES = { # expected structure
#     # Single reads with or without assembly
#     "ARGA00461": { # full pipeline
#         "reads": [f"{temp_input_longreads_folder}/ARGA00461/ARGA00461.fastq.gz"],
#         "assembly": [f"{temp_input_assembly_folder}/ARGA00461/ARGA00461_assembly.fasta"]
#     },
#     "ARGA00190": { # truncation cant be done
#         "reads": [f"{temp_input_longreads_folder}/ARGA00190/ARGA00190.fastq.gz"],
#         "assembly": []
#     },
#     # Paired reads with or without assembly
#     "ARGA00024": { # full pipeline
#         "reads": [f"{temp_input_shortreads_folder}/ARGA00024/ARGA00024_R1.fastq.gz", f"{temp_input_shortreads_folder}/ARGA00024/ARGA00024_R2.fastq.gz"],
#         "assembly": [f"{temp_input_assembly_folder}/ARGA00024/ARGA00024_assembly.fasta"]
#     },
#     "ARGA00025": { # truncation cant be done
#         "reads": [f"{temp_input_shortreads_folder}/ARGA00025/ARGA00025_R1.fastq.gz", f"{temp_input_shortreads_folder}/ARGA00025/ARGA00025_R2.fastq.gz"],
#         "assembly": [f"{temp_input_assembly_folder}/ARGA00025/ARGA00025_assembly.fasta"]
#     },
#     # Assembly without reads
#     "ARGA00031": { # everything with assembly
#         "reads": [],
#         "assembly": [f"{temp_input_assembly_folder}/ARGA00031/ARGA00031_assembly.fasta"]
#     },
# }