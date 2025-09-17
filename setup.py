# Python setup script to assist building/deploying the solver

import pathlib
import os
import typer
from rich import print

app = typer.Typer(help="Setup script to assist building/deploying the CNS solver")

@app.command()
def prepare_env() -> None:
    """ Prepares the environment for building and running the
    CNS solver"""

    # Ensure the environment variable CNS_SOLVER_HOME is set
    if "CNS_SOLVER_HOME" not in os.environ:
        current_dir = pathlib.Path(__file__).parent.resolve()
        os.environ["CNS_SOLVER_HOME"] = str(current_dir)
        print(f"Set CNS_SOLVER_HOME to {current_dir}")
    else:
        print(f"CNS_SOLVER_HOME is already set to {os.environ['CNS_SOLVER_HOME']}")

    # Source the .env file if it exists
    env_file = pathlib.Path(__file__).parent / ".env"
    if env_file.exists():
        os.system(f"source {env_file}")

    return

@app.command()
def build(build_type:str = "Release") -> int:
    """ Build the CNS solver """
    # Placeholder for build logic
    print("[blue]Building the CNS solver...[/blue]")

    build_dir = pathlib.Path(__file__).parent / "build"
    build_dir.mkdir(exist_ok=True)

    os.chdir(build_dir)
    istat = os.system(f"cmake -DCMAKE_BUILD_TYPE={build_type} ..")
    if istat != 0:
        print("CMake configuration failed.")
        return -1

    istat = os.system("make -j4")
    if istat != 0:
        print("[red]Build failed.[/red]")
        return -1
    os.chdir("..")
    print("[green]Build completed.[/green]")
    return 0

if __name__ == "__main__":
    app()
