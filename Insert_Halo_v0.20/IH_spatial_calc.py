import numpy as np
from IH_errors import StageException
from scipy.integrate import simpson

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
        """
        Returns x and y positions (in arcsec) for the selected halos. The spatial positions are generated depending on the spatial
        distribution setting defined by self.IHB.spatial which can be found in the global_parameters class found in Insert_halo_base
        This function requires that self.Nsub has already been generated.
        """
        dist_setting  = self.IHB.spatial

        if self.Nsub == None:
            raise StageException("The number of expected subhalos should be called before calling generate_spatial_distribution().")

        if dist_setting == "uniform":
            x, y = self.generate_uniform_distribution(self.Nsub, x_offset = self.IHB.center_x, y_offset = self.IHB.center_y)
            return x, y
        elif dist_setting == "einasto":
            x, y = self.generate_einasto_distribution(self.Nsub, x_offset = self.IHB.center_x, y_offset = self.IHB.center_y)
            return x, y 
        else:
            raise ValueError("Only a 'uniform' and 'einasto' spatial distribution settings are implemented so far")



    def generate_uniform_distribution(self, N, R = None, x_offset = 0, y_offset = 0):
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

        self.x_arcsec = x + x_offset
        self.y_arcsec = y + y_offset

        return self.x_arcsec, self.y_arcsec

    @staticmethod
    def einasto_prop(r, alpha = 0.678, r2 = 199):
        """
        Returns the proportionality Einasto profile which is proportional to n_sh(r) from https://arxiv.org/pdf/2403.16633
        Current version has alpha and r_{-2} testing values from the Aq-A-1 simulations of https://arxiv.org/pdf/0809.0898
        r2 is in units of kpc
        """
        return np.exp(-2 / alpha * (r / r2)** alpha)
    
    @staticmethod
    def normalize(r_array, function_values):
        A = simpson(function_values, r_array)
        return function_values / A
        
        
    def generate_einasto_distribution(self, N, R = None, x_offset = 0, y_offset = 0):
        """
        Generates spatial distribution positions (in arcsec) for N subhalos based on the Einasto profile
        This function draws values randomly using the Einasto profile to generate the radii values and then
        generates 3d positions before only taking the x and y positions and returning these
        """
        einasto_prop = spatial_distribution_generation.einasto_prop

        
        if R == None:
            R = self.IHB.cone_angle_arcsec / 2 * self.GM.kpc_per_arcsec()

        random_r_array = []

        for i in range(N):
            draw_found = False
            while not draw_found:
                drawr = np.random.rand()
                random_r = np.random.rand() * R
                comparison_p = einasto_prop(random_r)
                if comparison_p >= drawr:
                    random_r_array.append(random_r)
                    draw_found = True
                    
        random_r_array = np.array(random_r_array)/self.GM.kpc_per_arcsec()

        theta_random = np.random.uniform(0, np.pi, N)
        phi_random = np.random.uniform(0, 2 * np.pi, N)
        
        self.x_arcsec = random_r_array * np.sin(theta_random) * np.cos(phi_random) + x_offset 
        self.y_arcsec = random_r_array * np.sin(theta_random) * np.sin(phi_random) + y_offset
        return self.x_arcsec, self.y_arcsec
        

        
        







































        
        


