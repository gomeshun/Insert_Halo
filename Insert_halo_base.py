import numpy as np

class constants_and_conversion:
    """
    This class contains most relevant constants and means of converting their units. This class was heavily inspired by
    a similar class in Sashimi (https://github.com/shinichiroando/sashimi-c), however this class has slightly different definitions
    to work more easily with pyHalo instead.
    """
    def __init__(self):
        #Units and conversions:
        #Masses:
        self.Msun = 1
        self.gram = self.Msun/1.988e33
        self.kg = self.gram*1000

        #Distances:
        self.kpc = 1
        self.pc = self.kpc/1000
        self.AU = self.pc*np.tan(1/3600/180*np.pi)
        self.km = self.AU/(1.495978707e8)
        self.m = self.km/1000
        self.cm = self.m/100
        self.Mpc = self.kpc*1000

        #Times:
        self.s = 1
        self.day = self.s*3600*24
        self.year = self.day*365.242374
        self.Gyr = self.year*1e9

        #Angles:
        self.arcsec = 1
        self.arcmin = self.arcsec*60
        self.deg = self.arcmin*60
        self.rad = self.deg/np.pi * 180
        

        #Physical constants:
        
        #https://asd.gsfc.nasa.gov/Stephen.Merkowitz/G/Big_G.html
        self.G = 6.674215e-11 * self.m ** 3 / self.kg / self.s **2



class global_parameters(constants_and_conversion):
    """
    This class contains the parameters used globally by insert_halo. These include parameters that are 
    used by either pyHalo or Sashimi and/or used externally by insert_halo itself. Instantiating this 
    class requires the constants_and_conversion class, which contains physical constants and means of
    converting units and constants.
    """
    def __init__(self, zlens = 0.5, zsource = 2.0, DM_type = "CDM", sigma_0 = 20, vel_scale = 25, Mhost = 13, geometry_shape = "DOUBLE_CONE",
                shmfplindex = -1.9, msub_min = 6, msub_max = 10, sigma_sub = 0.025, cone_angle_arcsec = 6.0, host_halo_scale_factor = 0.88,
                redshift_scale_factor = 1.7, seed = 1234, randomization =  False):
        constants_and_conversion.__init__(self)
        self.zlens = zlens
        self.zsource = zsource
        self.DM_type = DM_type
        self.sigma_0 = sigma_0
        self.vel_scale = vel_scale
        self.Mhost = Mhost #Should be M_200 for Sashimi and should be in base log10 for pyHalo!


        # Relevant parameters for the geometry and subhalo mass function: 
        self.geometry_shape = geometry_shape
        self.shmfplindex = shmfplindex
        self.msub_min = msub_min
        self.msub_max = msub_max
        self.sigma_sub = sigma_sub
        self.cone_angle_arcsec = cone_angle_arcsec

        self.host_halo_scale_factor = host_halo_scale_factor # k1 in https://arxiv.org/pdf/1908.06983
        self.redshift_scale_factor = redshift_scale_factor # k2 in https://arxiv.org/pdf/1908.06983

        # Randomization settings:
        self.randomization = randomization
        self.seed = seed

        

class program_status:
    def __init__(self, current_status = "initializing"):
        self.current_status = current_status

    def status(self):
        """
        Returns current program status.
        """
        return self.current_status

    def set_status(self, new_status):
        """
        Sets new program status.
        """
        self.current_status = new_status
        
    


class Insert_halo_base(global_parameters):
    def __init__(self, zlens = 0.5, zsource = 2.0, DM_type = "CDM", sigma_0 = 20, vel_scale = 25, Mhost = 13, geometry_shape = "DOUBLE_CONE",
                shmfplindex = -1.9, msub_min = 6, msub_max = 10, sigma_sub = 0.025, cone_angle_arcsec = 6.0, host_halo_scale_factor = 0.88,
                redshift_scale_factor = 1.7, seed = 1234, randomization = False):
        global_parameters.__init__(self, zlens = zlens, zsource = zsource, DM_type = DM_type, sigma_0 = sigma_0,
                          vel_scale = vel_scale, Mhost = Mhost, geometry_shape = geometry_shape, shmfplindex = shmfplindex, 
                          msub_min = msub_min, msub_max = msub_max, sigma_sub = sigma_sub, cone_angle_arcsec = cone_angle_arcsec, 
                          host_halo_scale_factor = host_halo_scale_factor, redshift_scale_factor = redshift_scale_factor, seed = seed,
                          randomization = randomization)
        self.global_status = program_status()
        self.global_status.set_status("idle")





















