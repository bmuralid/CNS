"""
This requires invoke to be installed on the machine (not venv), which can be
done via:
    pipx install invoke
"""
import os
from datetime import datetime
from hashlib import md5
from pathlib import Path
from shutil import rmtree as shutil_rmtree
from typing import Optional

import invoke

PROJECT: str = "CNS"
SRC_PATH: Path = Path(__file__).parent
# VCPKG_TOOLCHAIN = SRC_PATH / "vcpkg/scripts/buildsystems/vcpkg.cmake"
WORKSPACE: Path = Path(__file__).parent
MD5: Optional[str] = None
BUILD_PATH: Optional[Path] = None
INSTALL_PATH: Optional[Path] = None
AMREX_ROOT_DIR = Path(__file__).parent/"./extern/amrex/install/nvhpc"

os.environ["AMREX_ROOT_DIR"] = str(AMREX_ROOT_DIR)

def get_md5(content: str) -> str:
    global MD5
    if MD5 is None:
        MD5 = md5(str.encode(content)).hexdigest()
    return MD5


def get_cmake_workspace() -> Path:
    # hash = get_md5(str(SRC_PATH))
    return WORKSPACE
    # return WORKSPACE / f"{PROJECT}_{SRC_PATH.name}_{hash}"


def get_build_path() -> Path:
    global BUILD_PATH
    if BUILD_PATH is None:
        BUILD_PATH = get_cmake_workspace() / "build"
    return BUILD_PATH


def get_install_path() -> Path:
    global INSTALL_PATH
    if INSTALL_PATH is None:
        INSTALL_PATH = get_cmake_workspace() / "install"
    return INSTALL_PATH


@invoke.task
def info(c, topic="all"):
    """Show project info."""
    if topic == "all":
        print(f"Project         = {PROJECT}")
        print(f"Source path     = {SRC_PATH}")
        print(f"Build path      = {get_build_path()}")
        print(f"Install path    = {get_install_path()}")
    elif topic == "build_path":
        print(get_build_path())
    elif topic == "install_path":
        print(get_install_path())
    else:
        print("Error: Valid 'topic' names are 'build_path'/'install_path'")

@invoke.task
def load_module(c):
    """ Load necessary modules."""
    module_list = [
        "mpi/latest"
    ]
    for mod in module_list:
        print(f"module load {mod}")


@invoke.task
def load_env(c):
    """Load environment variables from .env file if it exists."""
    env_file = SRC_PATH / ".env"
    if env_file.exists():
        # print(f"# Loading environment variables from {env_file}")
        with env_file.open("r") as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    key, value = line.strip().split("=", 1)
                    # Output in shell export format
                    print(f"export {key}={value}")
                    os.environ[key] = value
    else:
        print(f"# No .env file found at {env_file}, skipping.")

@invoke.task
def config(c, build_type:str="RelWithDebInfo", compiler:str="gcc-9"):
    """Run cmake configure."""
    do_config(c, build_type=build_type,compiler=compiler)


def do_config(c, build_type:str="RelWithDebInfo", compiler:str="gcc-9"):
    build_path = get_build_path()
    build_path.mkdir(parents=True, exist_ok=True)
    match compiler:
        case "gcc-9":
            # os.environ["CC"] = "gcc-9"
            # os.environ["CXX"] = "g++-9"
            # os.environ["FC"] = "gfortran-9"
            CC = "gcc-9"
            CXX = "g++-9"
            FC = "gfortran-9"

        case "nvhpc":
            CC = "nvc"
            CXX = "nvc++"
            FC = "nvfortran"
            os.environ["CC"] = "nvc"
            os.environ["CXX"] = "nvc++"
            os.environ["FC"] = "nvfortran"
        case _:
            CC = "gcc-9"
            CXX = "g++-9"
            FC = "gfortran-9"
            print(f"Warning: Unknown compiler '{compiler}'. Using system default.")
    cmd = [
        f"CXX={CXX}",
        f"CC={CC}",
        f"FC={FC}",
        "cmake",
        "-S",
        str(SRC_PATH),
        "-B",
        str(build_path),
        f"-DCMAKE_BUILD_TYPE={build_type}",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=1",
        "-DAMReX_GPU_BACKEND=CUDA",
    ]
    c.run(" ".join(cmd), pty=True)

    # Symlink compile_commands.json
    src_ccdb_file = SRC_PATH / "compile_commands.json"
    build_ccdb_file = build_path / "compile_commands.json"
    if build_ccdb_file.exists():
        if not src_ccdb_file.exists():
            src_ccdb_file.symlink_to(build_ccdb_file)


@invoke.task
def build(c, config=False):
    """Run builds via cmake."""
    build_path = get_build_path()

    if not build_path.exists():
        if config:
            do_config(c)
        else:
            print("Error: build path doesn't exist.")
            return

    cmd = ["cmake", "--build", str(build_path)]
    c.run(" ".join(cmd), pty=True)


@invoke.task
def install(c):
    """Run install via cmake."""
    build_path = get_build_path()
    install_path = get_install_path()

    if not build_path.exists():
        print("Error: build path doesn't exist.")
        return

    cmd = ["cmake", "--install", str(build_path)]
    c.run(" ".join(cmd), env={"DESTDIR": install_path}, pty=True)


@invoke.task
def clean(c):
    """Clean build directory."""
    # # Don't clean during CppCon and the week before/after
    # _, week, _ = datetime.now().isocalendar()
    # if week in (36, 37, 38):
    #     print(f"I'm sorry I can't do that Dave as the current week is {week}.")
    #     return

    build_path = get_build_path()
    if build_path.exists():
        shutil_rmtree(build_path)
        print(f"Cleaned {build_path}")
    else:
        print("Build path absent. Nothing to do.")


@invoke.task(pre=[clean])
def clean_all(c):
    """Clean build and install directory."""
    install_path = get_install_path()
    if install_path.exists():
        shutil_rmtree(install_path)
        print(f"Cleaned {install_path}")
    else:
        print("Install path absent. Nothing to do.")


@invoke.task
def ls(c):
    """List files using lsd and skip vcpkg and .cache folder"""
    cmd = [
        "lsd",
        "--tree",
        "--ignore-glob vcpkg",
        "--ignore-glob .cache",
    ]
    c.run(" ".join(cmd), pty=True)

@invoke.task
def ctags(c):
    """Generate ctags for the project."""
    cmd = [
        "ctags",
        "-R",
        "--languages=C,C++,Fortran",
        "--exclude=.git",
        "--exclude=build",
        "--exclude=extern",
        "-f",
        "tags",
        ".",
    ]
    c.run(" ".join(cmd), pty=True)

@invoke.task
def test(c):
    """Run tests via ctest."""
    build_path = get_build_path()

    if not build_path.exists():
        print("Error: build path doesn't exist.")
        return

    cmd = [
        "cd",
        f"{build_path} &&" ,
        "ctest",
        "--output-on-failure",
        "-C",
        "RelWithDebInfo",
        "-j",
        "4",
    ]
    c.run(" ".join(cmd), pty=True)
