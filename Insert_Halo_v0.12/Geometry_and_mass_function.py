import numpy as np
from scipy.integrate import simpson






def subhalo_mass_function(M_halo, z, m, sigma_sub, k1 = 0.88, k2 = 1.7, m0 = 1e8, alpha = -1.9):
    """
    Returns d^2N_sub/dm/dA as shown in https://arxiv.org/pdf/1908.06983, good plots to recreate: https://arxiv.org/pdf/1908.06983
    """

    log10F = k1 * np.log10(M_halo / 1e13) + k2 * np.log(z + 0.5)
    return sigma_sub/m0 * (m / m0) ** alpha * 10 ** log10F



def host_halo_scale(Mhalo, z, k1 = 0.88, k2 = 1.7):
    """
    Computes the host halo scaling factor F (https://arxiv.org/pdf/1908.06983) where k1 is the host scaling factor and k2 the redshift scaling factor
    """
    return 10**(k1 * np.log10(Mhalo / 1e13) + k2 * np.log10(z + 0.5))

def integrated_subhalo_mass_function(Mhalo, msub_min, msub_max, z, m, sigma_sub, k1 = 0.88, k2 = 1.7, m0 = 1e8, alpha = -1.9):
    """
    Returns N_sub / A from an integrated subhalo mass function as seen in https://arxiv.org/pdf/1908.06983. The returned value is the number of expected
    subhalos between the a mass of msub_min and msub_max at a redshift z. The A in the returned quantity is the area of the lens plane that we are considering.

    The area of the lens plane has to be evaluated externally from this function, but this has to be computed at the same redshift value as is given here.
    """
    int_factor = alpha + 1

    cfactor = sigma_sub * host_halo_scale(Mhalo, z, k1, k2) / (m0 ** int_factor)

    if int_factor != 0:
        integral = 1 / int_factor * (msub_max ** int_factor - msub_min ** int_factor)
    else:
        integral = np.log(msub_max / msub_min)

    return cfactor * integral


class geometry_and_massfunction():
    def __init__(self, SH_control_storage, cosmology):
        self.SH = SH_control_storage
        self.IHB = self.SH.IHB
        self.cosmology = cosmology
        

    def host_halo_scale(self):
        """
        Computes the host halo scaling factor F (https://arxiv.org/pdf/1908.06983) where k1 is the host scaling factor and k2 the redshift scaling factor
        """
        k1 = self.IHB.host_halo_scale_factor
        k2 = self.IHB.redshift_scale_factor
        return 10**(k1 * np.log10(10 ** self.IHB.Mhost / 1e13) + k2 * np.log10(self.IHB.zlens + 0.5))

    def integrated_subhalo_mass_function(self, m0 = 1e8):
        """
        Returns N_sub / A from an integrated subhalo mass function as seen in https://arxiv.org/pdf/1908.06983. The returned value is the number of expected
        subhalos between the a mass of msub_min and msub_max at a redshift z. The A in the returned quantity is the area of the lens plane that we are considering. The units of the returned value is [kpc^-2]
    
        The area of the lens plane has to be evaluated externally from this function, but this has to be computed at the same redshift value as is given here.
        """
        alpha = self.IHB.shmfplindex
        sigma_sub = self.IHB.sigma_sub
        Mhalo = self.IHB.Mhost
        msub_min = 10 ** self.IHB.msub_min # IHB.msub_min is in log10 so to get the actual mass we raise 10 to the power of IHB.msub_min
        msub_max = 10 ** self.IHB.msub_max # Same as for msub_mmin


        int_factor = alpha + 1
    
        cfactor = sigma_sub * self.host_halo_scale() / (m0 ** int_factor)
    
        if int_factor != 0:
            integral = 1 / int_factor * (msub_max ** int_factor - msub_min ** int_factor)
        else:
            integral = np.log(msub_max / msub_min)
    
        return cfactor * integral  

    def integrated_subhalo_mass_function_massbins(self, m_min, m_max, m0 = 1e8):
        """
        Returns N_sub / A from an integrated subhalo mass function as seen in https://arxiv.org/pdf/1908.06983. The returned value is the number of expected
        subhalos between the a mass of msub_min and msub_max at a redshift z. The A in the returned quantity is the area of the lens plane that we are considering. The units of the returned value is [kpc^-2]
    
        The area of the lens plane has to be evaluated externally from this function, but this has to be computed at the same redshift value as is given here.
        """
        alpha = self.IHB.shmfplindex
        sigma_sub = self.IHB.sigma_sub
        Mhalo = self.IHB.Mhost
        msub_min = 10 ** m_min # IHB.msub_min is in log10 so to get the actual mass we raise 10 to the power of IHB.msub_min
        msub_max = 10 ** m_max # Same as for msub_mmin


        int_factor = alpha + 1
    
        cfactor = sigma_sub * self.host_halo_scale() / (m0 ** int_factor)
    
        if int_factor != 0:
            integral = 1 / int_factor * (msub_max ** int_factor - msub_min ** int_factor)
        else:
            integral = np.log(msub_max / msub_min)
    
        return cfactor * integral  

    def lens_plane_area(self):
        """
        Computes the area of the lensplane at the redshift value defined by zlens and the cone opening angle of the geometry
        Furthermore this function uses the arcsec_per_kpc_proper from astropy using the defined cosmology to compute the 
        physical radius of the lensplane.

        To do: Check whether the geometry shape has any effect on the radius calculation in pyHalo (the pyHalo documentation
        makes it seem like this would be the case)
        """
        cone_angle = self.IHB.cone_angle_arcsec
        z = self.IHB.zlens
        #geometry_shape = self.IHB.geometry_shape
        LCDM_flat = self.cosmology.LCDM_flat

        """
        if geometry_shape == "DOUBLE_CONE":
            radius = 0.5 * cone_angle / LCDM_flat.arcsec_per_kpc_proper(z).value
        elif geometry_shape == "CYLINDER":
            # From pyHalo definition: 
            # the physical radius of the cylinder is the angular diameter distance to the lens times 2 * cone_opening_angle (in radians)
            # Cone angle in rad means we divide by self.IHB.rad 
            # the angular_diameter_distance function computes this in Mpc, we want this in kpc so we convert units
            radius = LCDM_flat.angular_diameter_distance(z).value * self.IHB.Mpc * cone_angle / self.IHB.rad
        elif geometry_shape == "CONE":
            radius = 0.5 * cone_angle / LCDM_flat.arcsec_per_kpc_proper(z).value
        else:
            raise ValueError(f"Defined geometry shape should be 'DOUBLE_CONE', 'CYLINDER' or 'CONE', but is {geometry_shape} in lens_plane_area function.")
        """
        radius = 0.5 * cone_angle / LCDM_flat.arcsec_per_kpc_proper(z).value

        # The computed radius is in kpc at the given redshift z = zlens so the Area is simply A = radius**2 * pi
        Area = radius**2 * np.pi
        return Area



    def kpc_per_arcsec(self, z = None, N = 1000):
        """
        Computes the number of kpc in one arcsec for a flat CDM cosmology
        """
    
        if z == None:
            z = self.IHB.zlens

        z_array = np.linspace(0, z, N)
        d_integrand = self.H0_H(z_array)
        d_integral = simpson(d_integrand, x = z_array)

        H0 = self.cosmology.H0
        c = self.IHB.c / self.IHB.km * self.IHB.s

        arcsec1 = 1 * self.IHB.arcsec
        rad1 = arcsec1 / self.IHB.rad

        D = rad1 * c / H0 * d_integral / (1 + z) * self.IHB.Mpc
        return D
        



    def H0_H(self, z):
        """
        Returns 1 / (H(z) / H0)
        """
        OmegaM = self.cosmology.OmegaM
        OmegaL = 1 - OmegaM 
        return 1 / np.sqrt(OmegaM * (1 + z)**3 + OmegaL)





















