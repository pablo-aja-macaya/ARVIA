
from pathlib import Path
import glob
from subprocess import run, PIPE
from pprint import pprint
from arvia.utils.console_log import CONSOLE_STDOUT, CONSOLE_STDERR, log_error_and_raise, rich_display_dataframe
import yaml
from colorama import Fore, Style
import traceback

def test_arvia_pipeline_input(main_output_folder: Path):
    # Variables
    # main_output_folder = "/home/usuario/Proyectos/Results/tests/arvia/test_pipeline"
    sratoolkit_url = "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.2.1/sratoolkit.3.2.1-ubuntu64.tar.gz"
    sratoolkit_targz = f"{main_output_folder}/sratoolkit.tar.gz"
    sratoolkit_folder = f"{main_output_folder}/sratoolkit"
    arvia_folder = f"{main_output_folder}/arvia"
    fasterq_dump = f"{main_output_folder}/sratoolkit/sratoolkit*/bin/fasterq-dump"
    input_yaml_file = f"{main_output_folder}/arvia_test_input.yaml"
    TEST_SAMPLES = {
        "SAMN44420630": {
            "SRR32735229": {
                "notes": "pacbio"
            },
            "SRR31104317": {
                "notes": "illumina"
            },
        }
    }

    CONSOLE_STDOUT.log(f"Running ARVIA input test..", style="info")

    # Assert folder does not exist
    try:
        assert not Path(main_output_folder).exists(), "Output folder exists, will not overwrite. Please delete or select non-existent folder and try again."
    except Exception as e:
        log_error_and_raise(f"{e}")

    # Create main directory
    CONSOLE_STDOUT.log("Creating output folder...")
    _ = run(f"mkdir -p {main_output_folder} {sratoolkit_folder} {arvia_folder}", check=True, shell=True)

    # Test ARVIA help menus
    try:
        CONSOLE_STDOUT.log("Testing ARVIA's help menus...")
        _ = run(f"arvia --help > {arvia_folder}/help.log", check=True, shell=True)
        _ = run(f"arvia run --help > {arvia_folder}/help.log", check=True, shell=True)
        _ = run(f"arvia test --help > {arvia_folder}/help.log", check=True, shell=True)

    except Exception as e:
        log_error_and_raise(f"{traceback.format_exc()}\nFailed running ARVIA's help messages: {e}")

    # Download sratoolkit and decompress
    try:
        CONSOLE_STDOUT.log("Downloading SRA Toolkit...")
        _ = run(f"wget -O {sratoolkit_targz} --quiet {sratoolkit_url} ", check=True, shell=True)
        CONSOLE_STDOUT.log("Decompressing SRA Toolkit...")
        _ = run(f"tar -xf {sratoolkit_targz} -C {sratoolkit_folder}", check=True, shell=True)
    except Exception as e:
        log_error_and_raise(f"{traceback.format_exc()}\nCould not download/decompress SRA Toolkit: {e}")

    # Download files
    d = {key: {"paired-end": [], "single-end": []} for key in TEST_SAMPLES.keys()}
    CONSOLE_STDOUT.log("Downloading files...")
    try:
        for biosample_id, v in TEST_SAMPLES.items():
            for sra_id, vv in v.items():
                CONSOLE_STDOUT.log(f"Downloading {sra_id} from {biosample_id}")
                sra_wd_folder = f"{main_output_folder}/dw/{biosample_id}/{sra_id}"
                _ = run(f"mkdir -p {sra_wd_folder}", check=True, shell=True)
                _ = run(f"{fasterq_dump} --quiet --split-3 --temp {main_output_folder} --outdir {sra_wd_folder} {sra_id} > {sra_wd_folder}/sra.log", check=True, shell=True)
                _ = run(f"pigz {sra_wd_folder}/*.fastq", check=True, shell=True)
                fs = glob.glob(f"{sra_wd_folder}/*.fastq.gz")
                if len(fs) not in [1,2]:
                    raise Exception(f"Unexpected. Number of files downloaded from {sra_id} ({biosample_id}) were not 1 or 2: {fs}")

                if len(fs) == 1:
                    d[biosample_id]["single-end"] = fs
                if len(fs) == 2:
                    d[biosample_id]["paired-end"] = fs
    except Exception as e:
        log_error_and_raise(f"{traceback.format_exc()}\nError downloading SRA file: {e}")


    # Generate input YAML for ARVIA
    CONSOLE_STDOUT.log("Generating input YAML...")
    with open(input_yaml_file, "w") as out_handle:
        for biosample_id, values in d.items():
            out_handle.write(f"# ---- {biosample_id} ----\n")  # section title

            yaml.dump(
                {
                    f"{biosample_id}_pe": {
                        "reads": values["paired-end"]
                    },
                }, 
                out_handle, default_flow_style=False, sort_keys=False
            )
            out_handle.write("\n")
            yaml.dump(
                {
                    f"{biosample_id}_se": {
                        "reads": values["single-end"]
                    },
                }, 
                out_handle, default_flow_style=False, sort_keys=False
            )

    # Test running ARVIA
    try:
        CONSOLE_STDOUT.log("Running ARVIA with input YAML and --previsualize...")
        _ = run(f"mkdir -p {arvia_folder}", check=True, shell=True)
        _ = run(f"arvia run --input_yaml {input_yaml_file} --output_folder {arvia_folder} --previsualize > {arvia_folder}/arvia.log", check=True, shell=True)

        # CONSOLE_STDOUT.log("Running ARVIA with input YAML...")
        # _ = run(f"arvia run --input_yaml {input_yaml_file} --output_folder {arvia_folder}", check=True, shell=True)

        CONSOLE_STDOUT.log("Test finished succesfully!", style="success")

    except Exception as e:
        log_error_and_raise(f"{traceback.format_exc()}\nFailed running ARVIA pipeline: {e}")

    return True

