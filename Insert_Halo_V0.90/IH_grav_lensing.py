import numpy as np
import copy
from scipy.interpolate import make_interp_spline

class Insert_Halo_Grav_Lensing():
    """
    This class contains the gravitational lensing computation code compatible with pyHalo and thus Insert Halo. All lensing computations
    possible with this class are done so through the gravitational lensing software Lenstronomy: https://github.com/lenstronomy/lenstronomy
    """
    def __init__(self):
        self.IHB = self.HC.IHB
        self.geometry = self.HC.GM
        if self.IHB.DM_type == "CDM":
            self.DM_structure = self.HC.PH.DM_dist
        elif self.IHB.DM_type == "SIDM":
            self.DM_structure = self.HC.PH.SIDM_dist
        #if self.HC.IHB.DM_type == "CDM":
        #    self.DM_structure = self.HC.PH.DM_dist
        #elif self.HC.IHB.DM_type == "SIDM":
        #    raise NotImplementedError("The SIDM model has not yet been implemented")
        self.DM_structure_adjusted = None
        import lenstronomy
        self.lenstronomy = lenstronomy
        from lenstronomy.LensModel.lens_model import LensModel
        self.LensModel = LensModel


    #####################################################
    """Lenstronomy and lensmodel functions"""
    #####################################################


    def prepare_lensmodel(self, mass_sheet_correction = True, include_subhalos = True, include_fieldhalos = True, include_host = True, host_lens_profile = "EPL", 
                         force_pointmasses = False):
        """
        Function to generate a Lenstronomy Lensmodel using the pyHalo realization data generated with Insert Halo. By default this generates a 
        Lensmodel using subhalos, field halos and the host halo, but these can each be turned off. The standard mass sheet correction from
        pyHalo can also be turned on or off here.

        ===============
        Parameters
        ===============
        mass_sheet_correction: Boolean to indicate whether or not to use the pyHalo mass sheet correction calculation (default: True)
        include_subhalos: Boolean to indicate whether to take the subhalo data into account (default: True)
        include_fieldhalos: Boolean to indicate whether to take the field halo data into account (default: True)
        include_host: Boolean to indicate whether to take the host halo data into account (default:True)
        """

        # This next condition causes this if statement to run only when include_subhalos is set to False or include_fieldhalos is set to False
        # But not when both are as this is taken care of later on in the code
        if (not include_subhalos or not include_fieldhalos) and not (not include_subhalos and not include_fieldhalos):
            keep_indices = []

            if not include_fieldhalos:
                for index, halo in enumerate(self.DM_structure.halos):
                    if halo.is_subhalo:
                        keep_indices.append(index)
            if not include_subhalos:
                for index, halo in enumerate(self.DM_structure.halos):
                    if not halo.is_subhalo:
                        keep_indices.append(index)

            
            # Making a deepcopy to only include the desired data while still keeping the original information in self.DM_structure
            self.DM_structure_adjusted = copy.deepcopy(self.DM_structure)
            self.DM_structure_adjusted.halos = np.array(self.DM_structure.halos)[keep_indices]
            self.DM_structure.halos = list(self.DM_structure.halos)
            self.DM_structure_adjusted.halos = list(self.DM_structure_adjusted.halos)

        
        elif include_subhalos and include_fieldhalos:
            self.DM_structure_adjusted = copy.deepcopy(self.DM_structure)

        if include_subhalos or include_fieldhalos:
            if not force_pointmasses:
                if self.IHB.DM_type == "CDM":
                    lens_model_list, lens_redshift_array, kwargs_subhalos, _ = self.DM_structure_adjusted.lensing_quantities(add_mass_sheet_correction = mass_sheet_correction)
                elif self.IHB.DM_type == "SIDM":
#self, density_profile_3d, kpc_per_arcsec, sigma_crit, gdecom_P = 8, ngaussian = 15, add_mass_sheet_correction = False, kwargs_mass_sheet = {}
                    lens_model_list, lens_redshift_array, kwargs_subhalos, _ = self.DM_structure_adjusted.lensing_quantities(density_profiles_3d = {"TNFWC_SH_tau": self.HC.TNFWC_SH_tau, "TNFWC_SH": self.HC.TNFWC_SH}, kpc_per_arcsec = self.geometry.kpc_per_arcsec_interp, sigma_crit = self.sigma_crit_lens, gdecom_P = self.IHB.gdecom_P, ngaussian = self.IHB.ngaussian, add_mass_sheet_correction = mass_sheet_correction, kwargs_mass_sheet = {})
            elif force_pointmasses:
                lens_model_list, lens_redshift_array, kwargs_subhalos = self.pointmasses_lensing()
                if mass_sheet_correction:
                    print("Warning, pointmass lensing does not account for mass sheet corrections!")
            #lensModel = LensModel(lens_model_list)
            if self.IHB.DM_type == "CDM":
                lensmodel_sub = self.LensModel(lens_model_list, self.IHB.zlens, self.IHB.zsource, lens_redshift_array, cosmo = self.DM_structure.astropy_instance, multi_plane = True)
            elif self.IHB.DM_type == "SIDM":
                lensmodel_sub = self.LensModel(lens_model_list, self.IHB.zlens, self.IHB.zsource, lens_redshift_array, cosmo = self.HC.PH.DM_dist.astropy_instance, multi_plane = True)                
        if include_host:
            if host_lens_profile == "EPL":
                self.host.profile = "EPL"
                theta_E, gamma_EPL = self.compute_EPL_from_NFW()
                host_lens_kwargs = {"theta_E": theta_E, "gamma": gamma_EPL, "e1": 0, "e2": 0, "center_x": self.IHB.center_x, "center_y": self.IHB.center_y}
                lensmodel_host = self.LensModel([self.host.profile])

                if include_subhalos or include_fieldhalos:
                    if self.IHB.DM_type == "CDM":
                        lensmodel_full = self.LensModel([self.host.profile] + lens_model_list, self.IHB.zlens, self.IHB.zsource, [self.IHB.zlens] + list(lens_redshift_array), cosmo = self.DM_structure.astropy_instance, multi_plane=True)
                    elif self.IHB.DM_type == "SIDM":
                        lensmodel_full = self.LensModel([self.host.profile] + lens_model_list, self.IHB.zlens, self.IHB.zsource, [self.IHB.zlens] + list(lens_redshift_array), cosmo = self.HC.PH.DM_dist.astropy_instance, multi_plane=True)                        
                    kwargs_all = [host_lens_kwargs] + kwargs_subhalos
                    return {"lm_host": lensmodel_host, "kw_host": [host_lens_kwargs], "lm_sub": lensmodel_sub, "kw_sub": kwargs_subhalos, "lm_all": lensmodel_full, "kw_all": kwargs_all}
                else:
                    return {"lm_host": lensmodel_host, "kw_host": [host_lens_kwargs]}
            
            elif host_lens_profile == "NFW":
                rs_angle, alpha_rs = self.compute_rs_alpha_rs()
                self.host.profile = "NFW"
                host_lens_kwargs = {"Rs": rs_angle, "alpha_Rs": alpha_rs, "center_x": self.IHB.center_x, "center_y": self.IHB.center_y}
                lensmodel_host = self.LensModel([self.host.profile])
                if include_subhalos or include_fieldhalos:
                    lensmodel_full = self.LensModel([self.host.profile] + lens_model_list, self.IHB.zlens, self.IHB.zsource, [self.IHB.zlens] + list(lens_redshift_array), cosmo = self.DM_structure.astropy_instance, multi_plane=True)
                    kwargs_all = [host_lens_kwargs] + kwargs_subhalos
                    return {"lm_host": lensmodel_host, "kw_host": [host_lens_kwargs], "lm_sub": lensmodel_sub, "kw_sub": kwargs_subhalos, "lm_all": lensmodel_full, "kw_all": kwargs_all}
                else:
                    return {"lm_host": lensmodel_host, "kw_host": [host_lens_kwargs]}   
            else:
                raise NameError(r"Only the EPL and NFW profile are currently implemented for host halo lensing, not {host_halo_profile} .")
            
        return {"lm_sub": lensmodel_sub, "kw_sub": kwargs_subhalos} 


#kwargs_macromodel = [{'theta_E': theta_E, 'center_x': host_x, 'center_y': host_y, 'e1': 0.1, 'e2': 0.05, 'gamma': 1.5},
#                    {'gamma1': 0.05, 'gamma2': 0.02}]
#lens_model_macro = ['EPL', 'SHEAR']


    def sigma_crit_lens(self, zlens = None, zsource = None):
        """
        Computes the critical lensing surface density value for lens and source redshifts zlens, zsource respectively

        sigma_crit = c^2 / (4 * pi * G) * D_s / (D_d * D_ds) 

        from a source like https://arxiv.org/pdf/2401.04165

        ===============
        Parameters
        ===============
        zlens: redshift value for the lens at which sigma_crit is to be computed, if the default value of None is used, this assumes the redshift value for
        the gravitational lens at zlens from the Insert Halo computation.
        zsource: redshift value for the source for which sigma_crit is to be computed. If the default value of None is used, this assume the redshift value for
        the source at zsource from the Insert Halo computation.
        """

        c = self.IHB.c / self.IHB.m * self.IHB.s # In units of m/s
        G = self.IHB.G / self.IHB.m**2 / self.IHB.kpc * self.IHB.Msun * self.IHB.s**2 # In units of m^2 kpc / (Msun * s^2)

        # Angular diameter distances are in units of kpc:

        if zlens == None and zsource == None:
            Ds = self.geometry.angular_diameter_distance(z_close = 0, z_far = self.IHB.zsource)
            Dd = self.geometry.angular_diameter_distance(z_close = 0, z_far = self.IHB.zlens)
            Dds = (self.geometry.angular_diameter_distance(z_close = 0, z_far = self.IHB.zsource, divide_one_z = False) - self.geometry.angular_diameter_distance(z_close = 0, z_far = self.IHB.zlens, divide_one_z = False)) / (1 + self.IHB.zsource)
        elif zlens != None and zsource != None:
            Ds = self.geometry.angular_diameter_distance(z_close = 0, z_far = zsource)
            Dd = self.geometry.angular_diameter_distance(z_close = 0, z_far = zlens)
            Dds = (self.geometry.angular_diameter_distance(z_close = 0, z_far = zsource, divide_one_z = False) - self.geometry.angular_diameter_distance(z_close = 0, z_far = zlens, divide_one_z = False)) / (1 + zsource)
        else:
            raise ValueError(f"Using the sigma_crit_lens function with manually given zlens and/or zsource requires that both zlens and zsource are given, currently zlens is {zlens} and zsource is {zsource}.")

        sigma_crit = c ** 2 / (4 * np.pi * G) * Ds / (Dd * Dds)

        #print(f"Ds: {Ds/1000:.2f}Mpc, Dd: {Dd/1000:.2f}Mpc, Dds: {Dds/1000:.2f}Mpc")
        
        return sigma_crit



    def compute_EPL_from_NFW(self):
        """
        This function approximates the parameter of an EPL profile (https://lenstronomy.readthedocs.io/en/latest/lenstronomy.LensModel.Profiles.html#module-lenstronomy.LensModel.Profiles.epl) for the host halo, assuming the host follows an NFW density profile. For this we first estimate the Einstein radius of
        the host using the equations found in https://arxiv.org/pdf/astro-ph/0112138v2 and then match the NFW convergence (https://arxiv.org/pdf/2310.03077)
        with the convergence for the EPL profile at the Einstein radius to find the EPL gamma parameter. This way we compute convergence slope of the NFW 
        profile matches that of the EPL profile at the Einstein radius, which for strong lensing is the most important.
        """
        Dd = self.geometry.angular_diameter_distance(z_close = 0, z_far = self.IHB.zlens)
        crit_lens = self.sigma_crit_lens()
        theta_arcsec = np.logspace(-3, 3, 1000) * self.IHB.arcsec
        theta_rad = theta_arcsec / self.IHB.rad
        r_theta = theta_rad * Dd
        #print(np.min(r_theta))
        #print(np.max(r_theta))


        host_properties = self.host.properties

        # Computing the average surface density for the given host halo:
        sigma_avg = self.avg_surface_density_NFW(r_theta, host_properties["rs"], host_properties["rhos"])

        # Numerically inverting the surface density fraction:
        x_f = make_interp_spline(np.flip(sigma_avg/crit_lens), np.flip(r_theta))

        # The Einstein radius can then be easily found by:
        theta_E = x_f(1) # In kpc
        
        kappa_NFW_theta_E = self.kappa_nfw(theta_E, host_properties["rs"], host_properties["rhos"])
        gamma_EPL = 3 - 2 * kappa_NFW_theta_E

        self.host.theta_E = theta_E / Dd * self.IHB.rad
        self.host.gamma_EPL = gamma_EPL

        #return theta_E / Dd * self.IHB.rad, gamma_EPL
        return self.host.theta_E, self.host.gamma_EPL
        
        
        
    def avg_surface_density_NFW(self, r, rs, rhos):
        """
        Function to compute the average surface density enclosed in r for an NFW profile, equation from https://arxiv.org/pdf/astro-ph/0112138v2
        """
        x = r/rs
        return 4 * rhos*rs * self.dl_surface_density_NFW(x) / x**2
    
    @staticmethod
    def dl_surface_density_NFW(x):
        """
        Computes g(x) from https://arxiv.org/pdf/astro-ph/0112138v2, the input x can be a single value or a numpy array.
        """
        if isinstance(x, np.ndarray):
            x = np.where((x < 1) + (x > 1), x, np.nan)
            xs1 = np.where(x < 1, x, np.nan)
            g_xs1 = np.where(np.isnan(xs1) ,np.nan , np.log(xs1 / 2) + 1 / np.sqrt(1 - xs1 **2) * np.arccosh(1 / xs1))
            xg1 = np.where(x > 1, x, np.nan)
            g_xg1 = np.where(np.isnan(xg1), g_xs1, np.log(xg1 / 2) + 1 / np.sqrt(xg1 ** 2 - 1) * np.arccos(1 / xg1))
    
            gxw1 = np.where(np.isnan(g_xs1), g_xg1, g_xs1)
            g_x = np.where(np.isnan(gxw1), 1 + np.log(1/2), gxw1)
            return g_x
    
        elif isinstance(x, int) or isinstance(x, float):
            if x < 1:
                return np.log(x / 2) + 1 / np.sqrt(1 - x **2) * np.arccosh(1 / x)
            elif x == 1:
                return 1 + np.log(1/2)
            elif x > 1:
                return np.log(x / 2) + 1 / np.sqrt(x ** 2 - 1) * np.arccos(1 / x)
            else:
                raise ValueError(r"This cannot be evaluated for an x value of {x}.")
        else:
            raise ValueError(r"The dl_avg_surface_density function only works for numpy arrays or single int or float values not {type(x)}")


    def kappa_nfw_dl(self, r, rs, corr = 1e-6):
        """
        Dimensionless convergence of the regular NFW profile, from https://iopscience.iop.org/article/10.1088/1538-3873/ac12db/pdf
        """
        u = r/rs
        
        # Adjusting the value of u where u == 1 as this causes devision by zero errors
        # Potentially to do: change this function such that for u == 1 F(x) returns 1/3
        u_corr = np.where(u == 1, 1 + corr, u)
    
        Fsmall1 = 1 / (np.sqrt(1 - u_corr **2)) * np.arctanh(np.sqrt(1 - u_corr ** 2))
        Flarge1 = 1 / (np.sqrt(u_corr **2 - 1)) * np.arctan(np.sqrt(u_corr **2 - 1))
        F = np.where(u_corr < 1, Fsmall1, Flarge1)
        return 1 / (u_corr **2 - 1) * (1 - F)
    
    
    def kappa_nfw(self, r, rs, rhos):
        """
        Convergence of the regular NFW density profile, from https://arxiv.org/pdf/2310.03077
        """
        return 2 * rhos*rs / self.sigma_crit_lens() * self.kappa_nfw_dl(r, rs)


    def kappa_EPL(self, theta, theta_E, gamma):
        """
        Returns the convergence of the EPL model (https://arxiv.org/pdf/1904.08400) given the einstein radius and power law slope.
        theta and theta_E should be in the same units (typically arcsec)
        """
        return (3 - gamma) / 2 * (theta_E / theta) ** (gamma - 1)



    def compute_rs_alpha_rs(self):
        """
        Computes the values of r_s and the deflection angle at r_s in units of [arcsec] for the use in the NFW host halo lensmodel of Lenstronomy.
        The equation to compute the deflection angle alpha_rs can be found in 
        https://lenstronomy.readthedocs.io/en/latest/lenstronomy.LensModel.Profiles.html#module-lenstronomy.LensModel.Profiles.nfw,
        however this equation is missing the critical lensing density sigma_crit, which we divide by in addition to also converting units.
        """
        host_props = self.host.properties
        Dd = self.geometry.angular_diameter_distance(z_close = 0, z_far = self.IHB.zlens)
        rs_arcsec = host_props["rs"] / Dd * self.IHB.rad

        crit_lens = self.sigma_crit_lens()
        alpha_rs = 4 * host_props["rs"] ** 2 * host_props["rhos"] * (1 + np.log(1/2)) / Dd / crit_lens * self.IHB.rad
        return rs_arcsec, alpha_rs


    def pointmasses_lensing(self):
        """
        Similar function to lensing_quantities from pyHalo, except here each halo (both subhalos and fieldhalos) are forced to be represented 
        using point masses.
        This function returns: lens_model_list, lens_redshift_array, kwargs_subhalos
        where lens_model_list is a list of Lenstronomy names for the type of model we are using (in this case 'PointMass')
        lens_redshift_array is a list of redshifts corresponding to each halo's redshift
        kwars_subhalos is a list of dictionaries for the properties of each halo (in this case, theta_E, center_x and center_y all in arcsec)

        Please note that mass correction is not implemented in this function.
        """

        lens_models = []
        redshifts = []
        kwargs_halos = []

        for halo in self.DM_structure_adjusted.halos:
            lens_models.append("POINT_MASS")
            redshifts.append(halo.z)
            if halo.is_subhalo:
                kwargs_halos.append({"theta_E": self.pointmass_einstein_radius(halo.z, halo.mass_evolved), "center_x": np.round(halo.x, 4), "center_y": np.round(halo.y, 4)})
            else:
                kwargs_halos.append({"theta_E": self.pointmass_einstein_radius(halo.z, halo.mass), "center_x": np.round(halo.x, 4), "center_y": np.round(halo.y, 4)})

        return lens_models, redshifts, kwargs_halos
        

        
    def pointmass_einstein_radius(self, halo_z, halo_mass):
        Ds = self.geometry.angular_diameter_distance(z_close = 0, z_far = self.IHB.zsource)
        Dd = self.geometry.angular_diameter_distance(z_close = 0, z_far = halo_z)
        Dds = (self.geometry.angular_diameter_distance(z_close = 0, z_far = self.IHB.zsource, divide_one_z = False) - self.geometry.angular_diameter_distance(z_close = 0, z_far = halo_z, divide_one_z = False)) / (1 + self.IHB.zsource)

        G = self.IHB.G / self.IHB.m**2 / self.IHB.kpc * self.IHB.s**2 * self.IHB.Msun
        c = self.IHB.c / self.IHB.m * self.IHB.s
        theta_E = np.sqrt(4 * G * halo_mass / c ** 2 * Dds / Dd / Ds)
        return theta_E * self.IHB.rad
        

    #####################################################
    """Subhalo convergence power spectrum functions"""
    #####################################################

    @staticmethod
    def rs_m_relation(mass, gamma = 1/3, m0 = 1e6, rs_0 = 0.1):
        """
        Equation (37) from https://link.aps.org/accepted/10.1103/PhysRevD.97.023001 to indicate the possible relation between r_s and the mass
        of subhalos.
        """
        rs = rs_0 * (mass / m0) ** gamma
        return rs















    