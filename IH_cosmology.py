import astropy.cosmology as astr_cosm
import numpy as np
import scipy.integrate as integrate


class IH_cosmology():
    """
    This class defines the revelant cosmological parameters (from Sashimi) and then uses those parameters to define a flat lambda CDM cosmology model
    using astropy.cosmology.FlatLambdaCDM
    """
    def __init__(self, SH_control_storage):
        self.SH = SH_control_storage
        self.IHB = self.SH.IHB

        self.H0 = 100*self.SH.sh_settings.h
        self.OmegaM = self.SH.sh_settings.OmegaM
        self.OmegaB = self.SH.sh_settings.OmegaB

        self.LCDM_flat = astr_cosm.FlatLambdaCDM(H0 = self.H0, Om0 = self.OmegaM, Ob0 = self.OmegaB)
        self.zlens_arcsec_per_kpc = self.LCDM_flat.arcsec_per_kpc_proper(self.IHB.zlens).value
        self.zlens_kpc_per_arcsec = 1/self.zlens_arcsec_per_kpc

    def get_virial_radius(self, mass ,z):
        """
        Computes the virial radius of a halo using equation (3) found in https://arxiv.org/pdf/1803.07691.
        First we find the multiplication factor Delta_c = mean_density / critical_density (https://arxiv.org/pdf/astro-ph/9710107)
        """
        d = self.OmegaM * (1 + z)**3 / (self.OmegaM * (1 + z)**3 + 1 - self.OmegaM) - 1
        #Delta_c = 18 * np.pi**2 + 82 * d - 39 * d**2
        Delta_c = 200
        rho_crit_z = self.SH.sh_settings.rhocrit(z) / self.SH.sh_settings.Msun * self.SH.sh_settings.kpc**3
        r_vir = (3 * mass / (4 * np.pi * Delta_c * rho_crit_z)) ** (1/3)
        return r_vir

    @staticmethod
    def SH_TNFW(r, rs, rt, rhos):
        x = r / rs
        dens = np.where(r <= rt, rhos / (x * (1 + x)**2), 0)
        return dens

    #def get_r200(self, rhos, rs, rt, z, rmax = 100, N = 1000):
    #    r_array = np.logspace(0, rmax*rs, N)
    #    rho_array = self.SH_TNFW(r_array, rs, rt, rhos)
    #    rho_crit_z = self.SH.sh_settings.rhocrit(z)
        
    def concentration_calc(self, mass , z):
        """
        Computes the concentration of a halo using the equation from the semi-analytical method described in https://arxiv.org/pdf/1502.00391 
        (Appendix B for a Plank cosmology)
        Mass should be in units of solar mass
        """

        if z <= 4:
            alpha = 1.7543 - 0.2766 * (1 + z) + 0.02039 * (1 + z) ** 2 
            beta = 0.2753 + 0.00351 * (1 + z) - 0.3038 * (1 + z) ** 0.0269 
            gamma = -0.01537 + 0.02102 * (1 + z) ** -0.1475
            return 10 ** (alpha + beta * np.log10(mass) * (1 + gamma * np.log10(mass) ** 2))
        elif z > 4:
            alpha = 1.3081 - 0.1078 * (1 + z) + 0.00398 * (1 + z) ** 2
            beta = 0.0223 - 0.0944 * (1 + z) ** -0.3907
            return 10 ** (alpha + beta * np.log10(mass))

    @staticmethod
    def calc_avg_density(density_profile, R_max, rs, rt, rhos, N = 1000):
        r_array = np.linspace(0.000001, R_max, N)
        rho_array = density_profile(r_array, rs, rt, rhos)
        rho_integrand = r_array ** 2 * rho_array
        rho_integral = integrate.simpson(y = rho_integrand, x = r_array)
        return 3 * rho_integral / R_max ** 3

    


            

#import astropy.cosmology as astr_cosm
#LCDM_flat = astr_cosm.FlatLambdaCDM(H0 = 70, Om0 = 0.3, Ob0 = 0.05)
#print(LCDM_flat.arcsec_per_kpc_proper(1))
#print(SH.integrated_subhalo_mass_function() * (SH.IHB.cone_angle_arcsec * 0.5 * 1/LCDM_flat.arcsec_per_kpc_proper(SH.IHB.zlens)).value**2 * np.pi)