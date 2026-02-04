import sys
from radprocess.interface.interface import launch_interface


def main():
    if sys.version_info < (3, 11):
        raise RuntimeError("radprocess requires Python ≥ 3.11")

    launch_interface()
