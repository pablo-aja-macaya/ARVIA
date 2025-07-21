import os
from colorama import Fore, Style
from arvia.utils.console_log import CONSOLE_STDOUT, CONSOLE_STDERR, log_error_and_raise

ARVIA_DIR = os.path.abspath(os.path.dirname(__file__))
WORKING_DIR = os.getcwd()
VERSION = "v0.1.0"

ascii =f"""
{Style.BRIGHT}{Fore.LIGHTYELLOW_EX}           _______      _______          {Style.RESET_ALL}
{Style.BRIGHT}{Fore.LIGHTYELLOW_EX}     /\   |  __ \ \    / /_   _|   /\    {Style.RESET_ALL}
{Style.BRIGHT}{Fore.LIGHTYELLOW_EX}    /  \  | |__) \ \  / /  | |    /  \   {Style.RESET_ALL}
{Style.BRIGHT}{Fore.LIGHTYELLOW_EX}   / /\ \ |  _  / \ \/ /   | |   / /\ \  {Style.RESET_ALL}
{Style.BRIGHT}{Fore.LIGHTYELLOW_EX}  / ____ \| | \ \  \  /   _| |_ / ____ \ {Style.RESET_ALL}
{Style.BRIGHT}{Fore.LIGHTYELLOW_EX} /_/    \_\_|  \_\  \/   |_____/_/    \_\{Style.RESET_ALL}
{Style.BRIGHT}{Fore.LIGHTYELLOW_EX}                                         {Style.RESET_ALL}
"""

print(ascii)


from arvia.utils.user_parser import get_parser
from arvia.utils.snakemake_common import run_snakemake


def main():
    if __name__ == "arvia.arvia":
        parser = get_parser()
        args = parser.parse_args()

        command = vars(args)["command"]
        parameters = vars(args)

        CONSOLE_STDOUT.log("Starting ARVIA...", style="info")
        if command == "run":
            if not parameters["reads"] and not parameters["assemblies"]:
                log_error_and_raise("Provide reads (--reads) and/or assemblies (--assemblies), please.")
            # Run
            run_snakemake(
                f"{ARVIA_DIR}/workflows/arvia.smk",
                parameters,
                "ARVIA",
            )

        # elif command == "db_install":
        #     run_snakemake(
        #         f"{ARVIA_DIR}/workflows/db_installation/db_installation.smk",
        #         parameters,
        #         "Modular database installation",
        #     )
