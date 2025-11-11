import numpy as np


class spatial_distribution_generation():
    def __init__(self, geometry_and_massfunction):
        self.GM = geometry_and_massfunction
        self.SH = self.GM.SH
        self.cosmology = self.GM.cosmology
        self.IHB = self.GM.IHB
        self.RNG = np.random.default_rng(self.IHB.seed)
        self.Nsub = None
        self.Nsub_random = None


    def calc_Nsub(self):
        """
        Computes the number of expected subhalos in the area at the main lens redshift according to the subhalo mass function
        """
        Nsub = self.GM.integrated_subhalo_mass_function() * self.GM.lens_plane_area()
        return int(round(Nsub, 0))


    def rand_Nsub(self):
        """
        Uses Poisson statistics to randomly draw a number given the expected number of subhalos
        """
        Nsub_ = self.calc_Nsub()
        self.Nsub = Nsub_
        random_Nsub = self.RNG.poisson(lam = Nsub_)
        self.Nsub_random = random_Nsub
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

        

        
        







































        
        


