import os

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
        

    def set_pymsesrc(self):
        """
        Write ~/.pymses/pymsesrc based on current configuration,
        then return the Pymsesrc object so Jupyter displays it.
        """
        ndust = ramses.read.hydro_file_descriptor(
            self.configparams.inout.ramses_output_dir
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