import numpy as np


class spatial_distribution_generation():
    def __init__(self, PH_run_storage):
        self.PH = PH_run_storage
        self.GM = self.PH.GM
        self.SH = self.GM.SH
        self.cosmology = self.GM.cosmology
        self.IHB = self.GM.IHB
        self.RNG = np.random.default_rng(self.IHB.seed)
        self.Nsub = None
        self.x_arcsec = None
        self.y_arcsec = None


    def calc_Nsub(self):
        """
        Computes the number of expected subhalos in the area at the main lens redshift according to the subhalo mass function
        """
        Nsub = self.GM.integrated_subhalo_mass_function() * self.GM.lens_plane_area()
        self.Nsub = int(round(Nsub, 0))
        return self.Nsub


    def rand_Nsub(self):
        """
        Uses Poisson statistics to randomly draw a number given the expected number of subhalos
        """
        Nsub_ = self.calc_Nsub()
        #self.Nsub = Nsub_
        random_Nsub = self.RNG.poisson(lam = Nsub_)
        self.Nsub = random_Nsub
        return random_Nsub

    def distribution_from_weights(self):
        if self.IHB.randomization:
            N_sub = self.rand_Nsub()
        else:
            N_sub = self.calc_Nsub()
        weights = self.SH.renormalize_weights(N_sub)
        probs = weights / N_sub

        dist_indices = self.RNG.choice(np.arange(len(probs)), size = N_sub, p = probs)

        return dist_indices

    def generate_spatial_distribution(self):
        dist_setting  = self.IHB.spatial

        if dist_setting == "uniform":
            pass

        else:
            raise ValueError("Only a 'uniform' spatial distributio setting is implemented so far")



    def generate_uniform_distribution(self, N, R = None):
        """
        Generates a uniform spatial distribution for N subhalos

        Maybe add the randomization seed here too for recreating results?
        """
        if R == None:
            R = self.IHB.cone_angle_arcsec / 2 
        
        theta = np.random.uniform(0, np.pi * 2, N)
        r = np.random.uniform(0, R ** 2, N)

        x = np.sqrt(r) * np.cos(theta)
        y = np.sqrt(r) * np.sin(theta)

        self.x_arcsec = x
        self.y_arcsec = y

        return x, y
        
        
        

        

        
        







































        
        


