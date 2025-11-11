import astropy.cosmology as astr_cosm


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





#import astropy.cosmology as astr_cosm
#LCDM_flat = astr_cosm.FlatLambdaCDM(H0 = 70, Om0 = 0.3, Ob0 = 0.05)
#print(LCDM_flat.arcsec_per_kpc_proper(1))
#print(SH.integrated_subhalo_mass_function() * (SH.IHB.cone_angle_arcsec * 0.5 * 1/LCDM_flat.arcsec_per_kpc_proper(SH.IHB.zlens)).value**2 * np.pi)