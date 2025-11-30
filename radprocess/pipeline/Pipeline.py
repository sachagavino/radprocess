import os
import json

from radprocess.pipeline.Grid import Grid
from radprocess.pipeline.Convert import Convert
from radprocess import radmc3d
from radprocess import ramses
from radprocess.utils.config import ConfigParams

from radprocess.constants.constants import au2cm, M_sun, R_sun

class Pipeline:

    def __init__(self):
        self.convert = Convert() 
        self.grid = Grid()
        self.configparams = ConfigParams()


    def read_pymsesrc(self):
        """
        Reads pymsesrc file from the current pymses.
        Returns a formatted string.
        """
        ramses_dir = self.configparams.ramsesoutput.ramses_output_dir
        try:
            nvar, variables, nb_dust = ramses.read.hydro_file_descriptor(ramses_dir)

            text = f"nvar = {nvar}\n\nVariables:\n"
            for idx, name in variables.items():
                text += f"  #{idx}: {name}\n"

            text += f"\nDust ratios detected: {nb_dust}\n"
            return text

        except Exception as e:
            return f"Error reading hydro_file_descriptor.txt:\n{e}"
        

    def read_hydro_descriptor(self):
        """
        Reads hydro_file_descriptor.txt from the current RAMSES directory.
        Returns a formatted string.
        """
        ramses_dir = self.configparams.ramsesoutput.ramses_output_dir
        try:
            nvar, variables, nb_dust = ramses.read.hydro_file_descriptor(ramses_dir)

            text = f"nvar = {nvar}\n\nVariables:\n"
            for idx, name in variables.items():
                text += f"  #{idx}: {name}\n"

            text += f"\nDust ratios detected: {nb_dust}\n"
            return text

        except Exception as e:
            return f"Error reading hydro_file_descriptor.txt:\n{e}"
        

    def read_sink_info(self):
        """
        Reads the sink_0000.info file from the RAMSES output directory
        using ramses.read.sink_info() function.
        Returns a string suitable for HTML display.
        """

        ramses_dir = self.configparams.ramsesoutput.ramses_output_dir
        try:
            sink_data = ramses.read.sink_info(ramses_dir)
        except Exception as e:
            return f"Error reading sink info: {e}"

        # Convert to aligned string
        columns = sink_data.columns
        rows = sink_data.rows

        lines = ["  ".join(columns)]
        for row in rows:
            lines.append("  ".join(str(row[col]) for col in columns))

        return "\n".join(lines)
        
    def read_pymsesrc(self):
        """Reads the ~/.pymses/pymsesrc JSON file and returns a dictionary."""
        pymses_directory = os.path.expanduser("~/.pymses")
        filename = os.path.join(pymses_directory, "pymsesrc")

        if not os.path.isfile(filename):
            raise FileNotFoundError(f"{filename} not found. Please create the pymsesrc file.")

        try:
            with open(filename, "r") as f:
                config_dict = json.load(f)
        except json.JSONDecodeError as e:
            raise RuntimeError(f"Error decoding JSON in {filename}:\n{e}")
        except Exception as e:
            raise RuntimeError(f"Error while reading {filename}:\n{e}")

        return config_dict

    def set_pymsesrc(self):
        """
        Write ~/.pymses/pymsesrc based on current configuration,
        then return the Pymsesrc object so Jupyter displays it.
        """
        ndust = ramses.read.hydro_file_descriptor(
            self.configparams.ramsesoutput.ramses_output_dir
        )[2]
        print(f'there is {ndust} dust species in the RAMSES simulation.')  
        # Write the file
        ramses.write.pymsesrc(self,
            ndust=ndust,
            rho=True,
            dustratios=self.configparams.pymsesrc.dustratios,
            vel=self.configparams.pymsesrc.vel,
            bl=self.configparams.pymsesrc.bl,
            br=self.configparams.pymsesrc.br,
            p=self.configparams.pymsesrc.p,
            xi=self.configparams.pymsesrc.xi,
            phi=self.configparams.pymsesrc.phi,
            g=self.configparams.pymsesrc.g,
        )

        # RETURN Pymsesrc config block → Jupyter displays it
        return self.configparams.pymsesrc
    
    def set_radmc3d_dir(self, path):
        """
        Ensure a RADMC-3D directory exists, store it in configparams,
        and optionally read wavelength/dust files if present.
        Returns a summary string for UI display.
        """
        # 1. Create directory if needed
        if not os.path.isdir(path):
            os.makedirs(path, exist_ok=True)
            created = True
        else:
            created = False

        # 2. Store the directory in config
        self.configparams.radoutput.radmc_output_dir = path

        messages = []
        if created:
            messages.append(f"Created RADMC-3D directory:\n  {self.configparams.radoutput.radmc_output_dir}")
        else:
            messages.append(f"Using existing RADMC-3D directory:\n  {self.configparams.radoutput.radmc_output_dir}")
    
        return "\n".join(messages)

        
    def thermal_radmc3d(self,
                        run=True, 
                        nphot=1e4, 
                        write_opac=True,
                        write_control=True, 
                        write_star=True, 
                        write_wave=True, 
                        write_mcmono=True, 
                        write_ext=True, 
                        **keywords):
        """ 
        Notes:
        run MC dust radiative transfer, open the resulting dust temperature as an array 
        and computes the surface-area weigthed temperature. If run == False, user assumes
        the RADMC3D output files already exist.
        -----
	    """	
        self.write_radmc3d(nphot_therm=nphot, 
                           write_opac=write_opac, 
                           write_control=write_control, 
                           write_star=write_star, 
                           write_wave=write_wave, 
                           write_mcmono=write_mcmono, 
                           write_ext=write_ext, 
                           **keywords)

        if write_control == False or write_opac == False or write_star == False or write_wave == False:
            print('WARNING: most RADMC3D input files will not be created. Will continue..\
                   but errors can be raised if one or more required input files are missing.\n')

        if run == True:
            self.run_thermal_radmc3d(nphot=nphot, **keywords)


    def write_radmc3d(self, 
                      write_dens, 
                      write_grid, 
                      write_opac, 
                      write_control, 
                      write_star, 
                      write_wave, 
                      write_mcmono, 
                      write_ext, 
                      **keywords):
        print('\nWRITING RADMC3D INPUT FILES:')
        print('----------------------------\n')
        #os.system("rm thermal/*.inp")
        if not os.path.exists('thermal'):
            os.makedirs('thermal')

        if write_control==True:
            radmc3d.write.control(**keywords)

    # LOAD RAMSES DATA
    def convert_rasmes2radmc(self):
        ramses_dir = self.configparams.ramsesoutput.ramses_output_dir
        radmc_dir = self.configparams.radoutput.radmc_output_dir

        clean = ramses_dir.strip().rstrip("/")   # remove whitespace + trailing slash
        folder = clean.rsplit("/", 1)[0] + "/"
        num = int(clean.rsplit("_", 1)[1])

        self.convert.to_radmc(folder, num, radmc_dir)