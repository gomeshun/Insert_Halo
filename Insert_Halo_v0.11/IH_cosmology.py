import astropy.cosmology as astr_cosm
import numpy as np


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
        Delta_c = 18 * np.pi**2 + 82 * d - 39 * d**2
        rho_crit_z = self.SH.sh_settings.rhocrit(z) / self.SH.sh_settings.Msun * self.SH.sh_settings.kpc**3
        r_vir = (3 * mass / (4 * np.pi * Delta_c * rho_crit_z)) ** (1/3)
        return r_vir





#import astropy.cosmology as astr_cosm
#LCDM_flat = astr_cosm.FlatLambdaCDM(H0 = 70, Om0 = 0.3, Ob0 = 0.05)
#print(LCDM_flat.arcsec_per_kpc_proper(1))
#print(SH.integrated_subhalo_mass_function() * (SH.IHB.cone_angle_arcsec * 0.5 * 1/LCDM_flat.arcsec_per_kpc_proper(SH.IHB.zlens)).value**2 * np.pi)