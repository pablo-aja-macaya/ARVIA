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


class CustomFormatter(
    RichHelpFormatter,
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.MetavarTypeHelpFormatter,
    argparse.RawTextHelpFormatter,
):
    pass


# Top-level parser
parser = argparse.ArgumentParser(
    allow_abbrev=False,
    formatter_class=ArgumentDefaultsRichHelpFormatter,
)
parser.add_argument(
    "--version",
    action="version",
    version=f"{Fore.YELLOW}ARVIA {VERSION}{Style.RESET_ALL}",
)

subparsers = parser.add_subparsers(help="command", required=True, dest="command")

#################################
# ---- Main ARVA ---- #
#################################
parser_run_arvia = subparsers.add_parser(
    "run",
    help=f"{Fore.YELLOW}ARVIA: Antibiotic Resistance Variant Identifier for Pseudomonas aeruginosa{Style.RESET_ALL}",
    description=textwrap.dedent(
        f"""
    {Fore.YELLOW}{Style.BRIGHT}ARVIA{Style.RESET_ALL}\n
    {Fore.WHITE}{Style.BRIGHT}Antibiotic Resistance Variant Identifier for Pseudomonas aeruginosa.{Style.RESET_ALL}
    """
    ),
    formatter_class=CustomFormatter,
)
parser_run_arvia__in_out = parser_run_arvia.add_argument_group("Input/Output")
parser_run_arvia__req_params = parser_run_arvia.add_argument_group("Required parameters")
parser_run_arvia__phylo = parser_run_arvia.add_argument_group("Phylogenetics")
parser_run_arvia__opt_params = parser_run_arvia.add_argument_group("Optional parameters")
# parser_run_arvia__conda_envs = parser_run_arvia.add_argument_group("Conda envs")

parser_run_arvia__in_out.add_argument(
    "--reads_folder",
    required=False,
    metavar="path",
    type=os.path.abspath,
    help=f"{Fore.LIGHTBLACK_EX}Input reads folder. (.fastq.gz/_R[1,2].fastq.gz/_[1,2].fastq.gz){Style.RESET_ALL}",
    dest="reads_folder",
)
parser_run_arvia__in_out.add_argument(
    "--assembly_folder",
    required=False,
    metavar="path",
    type=os.path.abspath,
    help=f"{Fore.LIGHTBLACK_EX}Input assembly folder (.fasta/.fna/.fa/.fas){Style.RESET_ALL}",
    dest="assembly_folder",
)
parser_run_arvia__in_out.add_argument(
    "--output_folder",
    required=True,
    metavar="path",
    type=os.path.abspath,
    help=f"{Fore.LIGHTBLACK_EX}Output folder{Style.RESET_ALL}",
    dest="output_folder",
)
parser_run_arvia__req_params.add_argument(
    "--previsualize",
    required=False,
    action="store_true",
    help=f"{Fore.LIGHTBLACK_EX}Previsualize pipeline to see if everything is correct{Style.RESET_ALL}",
    dest="previsualize",
)
parser_run_arvia__req_params.add_argument(
    "--cores",
    required=False,
    metavar="int",
    default=max(1, os.cpu_count()-1),
    type=int,
    help=f"{Fore.LIGHTBLACK_EX}Number of cores (default is available cores - 1){Style.RESET_ALL}",
    dest="cores",
)
parser_run_arvia__opt_params.add_argument(
    "--use_conda",
    required=False,
    action="store_true",
    help=f"{Fore.LIGHTBLACK_EX}If True, use conda environment specified by snakefile{Style.RESET_ALL}",
    dest="use_conda",
)
parser_run_arvia__opt_params.add_argument(
    "--barcodes",
    required=False,
    metavar="list",
    type=str,
    nargs="+",
    help=f"{Fore.LIGHTBLACK_EX}Space separated list of sample IDs. Only these samples will be processed{Style.RESET_ALL}",
    dest="barcodes",
)
parser_run_arvia__opt_params.add_argument(
    "--draw_wf",
    required=False,
    default=None,
    metavar="str",
    type=str,
    help=f"{Fore.LIGHTBLACK_EX}Draw pipeline to this path (PDF){Style.RESET_ALL}",
    dest="draw_wf",
)


# #################################
# # --- Database installation --- #
# #################################
# parser_db_installer = subparsers.add_parser(
#     "db_install",
#     help=f"{Fore.YELLOW}Modular database installation{Style.RESET_ALL}",
#     description=f"{Fore.YELLOW}{Style.BRIGHT}Install selected databases depending on user/pipeline needs{Style.RESET_ALL}",
#     formatter_class=CustomFormatter,
# )
# parser_db_installer.add_argument(
#     "--output_folder",
#     required=True,
#     metavar="path",
#     type=os.path.abspath,
#     help=f"{Fore.LIGHTBLACK_EX}Output folder{Style.RESET_ALL}",
#     dest="output_folder",
# )
# parser_db_installer.add_argument(
#     "--selected_pkgs",
#     required=True,
#     metavar="list",
#     type=str,
#     choices=[
#         "bactasys",
#         "bactasys-light",
#         "covpipe",
#         "metabactasys",
#         "all",
#     ],
#     nargs="+",
#     help=f"{Fore.LIGHTBLACK_EX}Pipeline packages (bactasys/covpipe/metabactasys/all){Style.RESET_ALL}",
#     dest="selected_db_packages",
# )
# parser_db_installer.add_argument(
#     "--previsualize",
#     required=True,
#     metavar="str",
#     type=str,
#     help=f"{Fore.LIGHTBLACK_EX}Previsualize pipeline to see if everything is correct{Style.RESET_ALL}",
#     choices=["yes", "no"],
#     dest="previsualize",
# )
# parser_db_installer.add_argument(
#     "--cores",
#     required=True,
#     metavar="int",
#     type=int,
#     help=f"{Fore.LIGHTBLACK_EX}Number of cores{Style.RESET_ALL}",
#     dest="cores",
# )
# parser_db_installer.add_argument(
#     "--use_conda",
#     required=False,
#     action="store_true",
#     help=f"{Fore.LIGHTBLACK_EX}If True, use conda environment specified by snakefile{Style.RESET_ALL}",
#     dest="use_conda",
# )
# parser_db_installer.add_argument(
#     "--draw_wf",
#     required=False,
#     default=None,
#     metavar="str",
#     type=str,
#     help=f"{Fore.LIGHTBLACK_EX}Draw pipeline to this path (PDF){Style.RESET_ALL}",
#     dest="draw_wf",
# )



def get_parser(parser=parser, subparsers=subparsers):
    return parser
