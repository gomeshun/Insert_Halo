import numpy as np
from scipy.interpolate import make_interp_spline
from scipy.optimize import least_squares
import math
from scipy.integrate import simpson


class Insert_Halo_SIDM_Distribution():
    """
    Class to hold the dark matter distribution as opposed to pyHalo. This class is specifically meant to be used with the SIDM part of Insert_Halo. For CDM
    Insert_Halo uses the pyHalo dark matter realization class.
    """
    def __init__(self):
        self.halo_list = []

    def halos(self):
        return self.halo_list

    def subhalos(self):
        subhalos = []
        for halo in self.halo_list:
            if halo.is_subhalo:
                subhalos.append(halo)
        return subhalos

    def field_halos(self):
        field_halos = []
        for halo in self.halo_list:
            if not halo.is_subhalos:
                field_halos.append(halo)
        return field_halos

    def lensing_quantities(self, density_profile_3d, kpc_per_arcsec, sigma_crit, gdecom_P = 8, ngaussian = 15, add_mass_sheet_correction = False, kwargs_mass_sheet = {}):
        """
        Lensing_quantities function as an alternative to the one in pyHalo of the same name. Here we convert the SIDM halos from Sashimi into Gaussians using the 
        Gaussian decomposition algorithm. Should return:
        
        lens_model_list, redshift_array, kwargs_lens, numerical_interp
        """

        lens_model_list = []
        redshift_array = []
        kwargs_lens_list = []
        numerical_interp = None  # This is unused in Insert_Halo, however it is still used here as an easy way to generate a lensmodel in Lenstronomy later without having
        # to distinguish between this particular class and the equivalent pyHalo one.
        GaussD = gaussian_decomposition(gdecom_P, ngaussian)
        

        for halo in self.halo_list:
            lens_model_list.append(["MULTI_GAUSSIAN"])
            redshift_array.append([halo.z])
            profile_params = halo.profile_args
            r_array, sigma = self.lens_projected_density(density_profile_3d, profile_params, np.log10(10 * profile_params["rt"]))
            proj_dens_interp = make_interp_spline(r_array, sigma)

            analytical_interp, analytical_interp_params = self.fit_projected_interp(r_array, proj_dens_interp(r_array), profile_params["rc"], profile_params["rs"]
                                                                                   ,profile_params["rt"])
            logsigma_min = np.log10(1e-3 * profile_params["rc"])
            logsigma_max = np.log10(10 * profile_params["rt"])
            gaussian_amplitudes, gaussian_sigmas, _, = GaussD.decompose(logsigma_min, logsigma_max, analytical_interp)

            pos_cond = gaussian_amplitudes > 0
            gaussian_sigmas_arcsec = gaussian_sigmas / kpc_per_arcsec
            gaussian_amplitudes_kappa = gaussian_amplitudes / sigma_crit * 2* np.pi * gaussian_sigmas_arcsec**2

            kwargs_lens_list.append({"amp": list(gaussian_amplitudes[pos_cond]), "sigma": list(gaussian_sigmas[pos_cond]), "center_x": halo.x, "center_y": halo.y, "scale_factor": 1})

        return lens_model_list, redshift_array, kwargs_lens_list, numerical_interp


    
    def lens_projected_density(self, density_profile_3d, density_profile_params, log_rmax, Nr = 1000, los_max = 1e3, losN = 1000):
        r_array = np.logspace(-3, log_rmax, Nr)

        los_array = np.linspace(-los_max, los_max, losN)
        los_array_expanded = np.expand_dims(los_array, axis = -1)
        los_array_2d = np.repeat(los_array_expanded, Nr, axis = 1)
        dplos_array = np.sqrt(los_array_2d ** 2 + r_array **2)
        integrand = density_profile_3d(dplos_array, density_profile_params)
        integral = simpson(integrand, x = los_array_2d, axis = 0)
        return r_array, integral
    
    def projected_interp_analytical(self, radius, amplitude, rc, rs, rt):
        radius = np.asarray(radius)
        term_core = np.sqrt(radius ** 2 + rc ** 2)
        term_scale = np.sqrt(radius ** 2 + rs ** 2)
        term_trunc = rt ** 2 / (radius ** 2 + rt ** 2)
        return amplitude / (term_core * term_scale) * term_trunc
    
    
    def fit_projected_interp(self, r_values, profile_values, rc_guess, rs_guess, rt_guess):
        r_values = np.asarray(r_values, dtype=float)
        profile_values = np.clip(np.asarray(profile_values, dtype=float), 1e-300, None)
        initial = np.array([
            np.log(np.max(profile_values) * np.sqrt(r_values[0] ** 2 + rc_guess ** 2) * np.sqrt(r_values[0] ** 2 + rs_guess ** 2)),
            np.log(rc_guess),
            np.log(rs_guess),
            np.log(rt_guess),
        ])
    
        def residuals(theta):
            amplitude, rc, rs, rt = np.exp(theta)
            model = self.projected_interp_analytical(r_values, amplitude, rc, rs, rt)
            return np.log(model) - np.log(profile_values)
    
        fit_result = least_squares(residuals, initial, max_nfev=4000, xtol=1e-12, ftol=1e-12, gtol=1e-12)
        amplitude, rc, rs, rt = np.exp(fit_result.x)
    
        def evaluate(radius):
            return self.projected_interp_analytical(radius, amplitude, rc, rs, rt)
    
        return evaluate, {
            "amplitude": amplitude,
            "rc": rc,
            "rs": rs,
            "rt": rt,
            "cost": fit_result.cost,
        }

    def Gaussian(self, r, amplitude, sigma):
        return amplitude * np.exp(-r**2 / (2 * sigma**2))
    
    def reconstruct_profile(self, r, amplitudes, sigmas):
        sum_of_gauss = 0
        for a, sigma in zip(amplitudes, sigmas):
            if a >= 0:
                sum_of_gauss += self.Gaussian(r, a, sigma)
    
        return sum_of_gauss






class gaussian_decomposition():
    """
    This class is meant to decompose an arbitrary halo density profile according to the method described in https://arxiv.org/pdf/1906.08263 into Gaussians
    for the purpose of being able to use it for strong gravitational lensing analysis.
    """
    def __init__(self, P, NG):
        self.P = P
        self.NG = NG
        self.xi_values = np.zeros(2*self.P + 1, dtype = np.float64)
        self.weights = np.zeros(2*self.P + 1, dtype = np.float64)        

        self.node_values = np.zeros(2*self.P + 1, dtype = np.complex64)
        self.precompute_nodes()
        self.precompute_weights()



    def decompose(self, logsigma_min, logsigma_max, num_interp_function):
        """
        Actual fuction for the Gaussian decomposition, logsigma_min is the smallest standard deviation value considered and logsigma_max the maxium standard deviation.
        num_interp_function is the numerical interpolation function of the projected density of the halo density profile.
        """
        sigma_array = np.logspace(logsigma_min, logsigma_max, self.NG)

        weights_expanded = np.expand_dims(self.weights, axis = -1)
        weights_2d = np.repeat(weights_expanded, self.NG, axis = 1)

        nodes_expanded = np.expand_dims(self.node_values, axis = -1)
        nodes_2d = np.repeat(nodes_expanded, self.NG, axis = 1)

        nodes_sigma = nodes_2d * sigma_array

        f_sigma_ = weights_2d * np.real(num_interp_function(nodes_sigma))
        #profile_i = self.profile(nodes_sigma)

        f_sigma = np.sum(f_sigma_, axis = 0)
        sigma_lin = np.linspace(logsigma_min, logsigma_max, self.NG)
        sigma_lin_spacing = np.mean(sigma_lin[1:] - sigma_lin[:-1])
        

        
        weights_int = np.ones(self.NG)
        weights_int[0] = 1/2
        weights_int[-1] = 1/2

        # 10**(self.P / 3) is from the constant term in the weights, which we add here and 1 / log10(e) is from converting the logarithmically
        # spaced sigma values to natural log from log10
        amplitudes = weights_int * f_sigma * sigma_lin_spacing / np.sqrt(2 * np.pi) * 10**(self.P / 3) / np.log10(np.exp(1)) 
        return amplitudes, sigma_array, f_sigma

    def precompute_nodes(self):
        for n in range(2*self.P + 1):
            self.node_values[n] = self.node(n)
    def node(self, n):
        first_term = 2 * self.P * np.log(10) / 3
        second_term = 2 * np.pi * np.complex64(1.j) * n
        return np.sqrt(first_term + second_term)


    def precompute_weights(self):
        for n in range(2*self.P + 1):
            self.xi_values[-(n + 1)] = self.weight(2 * self.P - n)
            #self.weights[-(n + 1)] = (-1)**(self.P * 2 - n) * 2*np.sqrt(2*np.pi)*10**(self.P / 3)*self.xi_values[-(n + 1)]
            self.weights[-(n + 1)] = (-1)**(self.P * 2 - n) * 2*np.sqrt(2*np.pi)*self.xi_values[-(n + 1)]

    def weight(self, n):
        if n == 2*self.P:
            return 2**(-self.P)
        elif n < self.P*2 and n > self.P:
            return self.xi_values[-(2*self.P - n)] + 2**(-self.P) * math.comb(self.P, 2*self.P - n)
        elif n <= self.P and n >= 1:
            return 1
        elif n == 0:
            return 1/2

            

class SIDM_halo():
    """
    Alternative SIDM halo class. This is a class which (potentially) replaces the SIDM halos implemented in pyHalo. The reason for this is that 
    the SIDM halo density profiles used in Insert_Halo are decomposed into Gaussian for the gravitational lensing analysis. Since this is not 
    compatible with pyHalo this alternative class is implemented here.
    """
    def __init__(self, mass, x, y, z, sub_flag, profile_args, has_collapsed, t_over_tc):
        self.mass = mass
        self.x = x
        self.y = y
        self.z = z
        self.is_subhalo = sub_flag
        self.profile_args = profile_args
        self.has_collapsed = has_collapsed
        self.t_over_tc = t_over_tc