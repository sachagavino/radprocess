import sys
import numpy as np
from radprocess.interface.interface import launch_interface

if __name__ == "__main__":

    if sys.version_info < (3, 11):
        raise RuntimeError("radprocess requires Python ≥ 3.11")

    if np.__version__ != "1.26.4":
        raise RuntimeError(
            f"Unsupported NumPy version {np.__version__}. "
            "Use numpy==1.26.4"
        )
        
    launch_interface()