import sys, os, inspect

from . Pipeline import Pipeline

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 


class Interface(Pipeline):
    def __init__(self, path_to_ramses_model, *args, **kwargs):
        """
        Initialize the Interface class.
        Parameters:
        -----------
        path_to_ramses_model : str
            Path to the RAMSES model.
        *args, **kwargs : Additional arguments for the Pipeline class.
        """
        # Initialize the parent class (Pipeline)
        super().__init__(*args, **kwargs)

        # Store the path to the RAMSES model
        self.ramses_path = path_to_ramses_model

    def set_pymsesrc(self, dustratios=True, rho=True, vel=True, Br=False):
        """
        Set the pymsesrc file for the RAMSES data.
        Parameters:
        -----------
        dustratios : bool
            Whether to include the dust ratios field.
        rho : bool
            Whether to include the density field.
        vel : bool
            Whether to include the velocity field.
        Br : bool
            Whether to include the radial magnetic field component.
        """
        self.pymsesrc = self.convert.update_pymsesrc(dustratios=dustratios, rho=rho, vel=vel, Br=Br)

    def do_ramses2radmc(self):
        self.grid.add_radmc_grid(self.convert.to_radmc(self.ramses_path))


    def do_ramses2polaris(self):
        self.grid.add_polaris_grid(self.convert.to_polaris(self.ramses_path))

    def do_radmc2polaris(self):
        self.grid.add_convert_grid(self.convert.radmc_to_polaris(self.ramses_path))