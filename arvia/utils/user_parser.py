from colorama import Fore, Style
from pathlib import Path
import argparse, textwrap
import os
import git

# import arvia
from arvia.arvia import VERSION
from rich_argparse import RichHelpFormatter, ArgumentDefaultsRichHelpFormatter
from arvia.utils.console_log import CONSOLE_STDOUT, CONSOLE_STDERR


RichHelpFormatter.styles["argparse.groups"] = "white bold"
RichHelpFormatter.styles["argparse.default"] = ""
# RichHelpFormatter._max_help_position = 52
# RichHelpFormatter._width = 200

class CustomFormatter(
    RichHelpFormatter,
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.MetavarTypeHelpFormatter,
    argparse.RawTextHelpFormatter,
):
    pass


# Top-level parser
parser = argparse.ArgumentParser(
    description=f"{Fore.YELLOW}{Style.BRIGHT}ARVIA:{Style.RESET_ALL} Antibiotic Resistance Variant Identifier for Pseudomonas aeruginosa",
    allow_abbrev=False,
    formatter_class=RichHelpFormatter,
    usage=argparse.SUPPRESS
)
parser.add_argument(
    "-v", "--version",
    action="version",
    version=f"ARVIA {VERSION}",
)

subparsers = parser.add_subparsers(help="command", required=True, dest="command")

########################
# ---- Main ARVIA ---- #
########################
parser_run_arvia = subparsers.add_parser(
    "run",
    help=f"Run ARVIA",
    description=f"{Fore.YELLOW}{Style.BRIGHT}ARVIA:{Style.RESET_ALL} Antibiotic Resistance Variant Identifier for Pseudomonas aeruginosa",
    formatter_class=lambda prog: CustomFormatter(prog, max_help_position=60, width=140),
)
parser_run_arvia__in_out = parser_run_arvia.add_argument_group("Input/Output")
parser_run_arvia__req_params = parser_run_arvia.add_argument_group("Required parameters")
parser_run_arvia__opt_params = parser_run_arvia.add_argument_group("Optional parameters")

parser_run_arvia__in_out.add_argument(
    "-r", "--reads", 
    required=False,
    metavar="path",
    type=os.path.abspath,
    nargs="+",
    help=f"Input reads files. Can be paired-end or single-end and must follow one of these structures: '{{sample_id}}.fastq.gz' / '{{sample_id}}_R[1,2].fastq.gz' / '{{sample_id}}_[1,2].fastq.gz' / '{{sample_id}}_S\\d+_L\\d+_R[1,2]_\\d+.fastq.gz'",
    dest="reads",
)
parser_run_arvia__in_out.add_argument(
    "-a", "--assemblies", 
    required=False,
    metavar="path",
    type=os.path.abspath,
    nargs="+",
    help=f"Input assembly files. Must follow one of these structures: '{{sample_id}}.{{fasta,fna,fa,fas}}'",
    dest="assemblies",
)
parser_run_arvia__in_out.add_argument(
    "-o", "--output_folder", 
    required=True,
    metavar="path",
    type=os.path.abspath,
    help=f"Output folder",
    dest="output_folder",
)
parser_run_arvia__opt_params.add_argument(
    "-c", "--cores", 
    required=False,
    metavar="int",
    default=max(1, os.cpu_count()-1),
    type=int,
    help=f"Number of cores (default is available cores - 1)",
    dest="cores",
)
parser_run_arvia__opt_params.add_argument(
    "-p", "--previsualize", 
    required=False,
    action="store_true",
    help=f"Previsualize pipeline to see if everything is as expected",
    dest="previsualize",
)
parser_run_arvia__opt_params.add_argument(
    "--use_conda",
    required=False,
    action="store_true",
    help=f"If True, use conda environment specified by snakefile",
    dest="use_conda",
)
parser_run_arvia__opt_params.add_argument(
    "--barcodes",
    required=False,
    metavar="str",
    type=str,
    nargs="+",
    help=f"Space separated list of sample IDs. Only these samples will be processed",
    dest="barcodes",
)
parser_run_arvia__opt_params.add_argument(
    "--draw_wf",
    required=False,
    default=None,
    metavar="str",
    type=str,
    help=f"Draw pipeline to this path (PDF",
    dest="draw_wf",
)

########################
# ---- Test ARVIA ---- #
########################
parser_test_arvia = subparsers.add_parser(
    "test",
    help=f"Test tool with a set of test files to ensure it is working properly.",
    description=f"Test {Fore.YELLOW}{Style.BRIGHT}ARVIA{Style.RESET_ALL} with a set of test files to ensure it is working properly",
    formatter_class=lambda prog: CustomFormatter(prog, max_help_position=60, width=140),
)

parser_test_arvia__in_out = parser_test_arvia.add_argument_group("Input/Output")
parser_test_arvia__req_params = parser_test_arvia.add_argument_group("Required parameters")
parser_test_arvia__opt_params = parser_test_arvia.add_argument_group("Optional parameters")

parser_test_arvia__in_out.add_argument(
    "-o", "--output_folder", 
    required=True,
    metavar="path",
    type=os.path.abspath,
    help=f"Output folder",
    dest="output_folder",
)
parser_test_arvia__opt_params.add_argument(
    "-c", "--cores", 
    required=False,
    metavar="int",
    default=max(1, os.cpu_count()-1),
    type=int,
    help=f"Number of cores (default is available cores - 1)",
    dest="cores",
)
parser_test_arvia__opt_params.add_argument(
    "-p", "--previsualize", 
    required=False,
    action="store_true",
    help=f"Previsualize pipeline to see if everything is as expected",
    dest="previsualize",
)
parser_test_arvia__opt_params.add_argument(
    "--use_conda",
    required=False,
    action="store_true",
    help=f"If True, use conda environment specified by snakefile",
    dest="use_conda",
)
parser_test_arvia__opt_params.add_argument(
    "--barcodes",
    required=False,
    metavar="str",
    type=str,
    nargs="+",
    help=f"Space separated list of sample IDs. Only these samples will be processed",
    dest="barcodes",
)
parser_test_arvia__opt_params.add_argument(
    "--draw_wf",
    required=False,
    default=None,
    metavar="str",
    type=str,
    help=f"Draw pipeline to this path (PDF",
    dest="draw_wf",
)

def get_parser(parser=parser, subparsers=subparsers):
    return parser
