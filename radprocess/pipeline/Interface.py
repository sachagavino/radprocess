import sys
import os
import inspect

from radprocess.pipeline.Pipeline import Pipeline
#from radprocess import pymses
from radprocess import radmc3d
from radprocess import ramses
#from radprocess.utils.config import ConfigParams
from radprocess.ramses.read import sink_info

# currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
# parentdir = os.path.dirname(currentdir)
# sys.path.insert(0,parentdir) 


class Interface(Pipeline):
    # def __init__(self, *args, **kwargs):
    #     """
    #     Initialize the Interface class.
    #     Parameters:
    #     -----------
    #     path_to_ramses_model : str
    #         Path to the RAMSES model.
    #     *args, **kwargs : Additional arguments for the Pipeline class.
    #     """
    #     # Initialize the parent class (Pipeline)
    #     super().__init__(*args, **kwargs)

    #     #self.configparams = ConfigParams()
    #     #self.sinkinfo = sink_info(self.configparams.inout.)
    #     #self.io_params = self.configparams.inout
    #     #self.sim_params = self.configparams.sim

    @property
    def sinkinfo(self):
        """
        Load the sink_XXXXX.info file from the RAMSES output directory.
        """

        # 1) Get directory from config
        path = self.configparams.inout.ramses_output_dir

        # Ensure trailing slash
        if not path.endswith("/"):
            path += "/"

        # Check RAMSES directory
        if not os.path.isdir(path):
            raise FileNotFoundError(
                f"RAMSES directory not found:\n    {path}\n"
                f"Please update: model.configparams.inout.ramses_output_dir"
            )
        return sink_info(path)
        #self.sinkinfo = sink_info(self.configparams.inout.ramses_output_dir)  

    def do_ramses2radmc(self):
        self.grid.add_radmc_grid(self.convert.to_radmc(self.ramses_path))

    def do_ramses2polaris(self):
        self.grid.add_polaris_grid(self.convert.to_polaris(self.ramses_path))

    def do_radmc2polaris(self):
        self.grid.add_convert_grid(self.convert.radmc_to_polaris(self.ramses_path))