from scipy.interpolate import make_interp_spline
import numpy as np
from scipy.integrate import simpson
from IH_errors import StageException

class SIDM_computations():
    def __init__(self, geometry_and_massfunction):
        self.GM = geometry_and_massfunction
        self.sh = self.GM.SH.sh_settings
        self.IHB = self.GM.IHB

        self.sigma_eff_interp = self.sh.sigma_eff_m_interpolate_analytical(sigma0_m = self.IHB.sigma_0 * self.sh.cm **2 / self.sh.gram,
                                                                           w = self.IHB.vel_scale * self.sh.km / self.sh.s)
        self.collapse_time = self.sh.t_collapse


#sigma_eff_interp = IHf.base_run.sh.sigma_eff_m_interpolate_analytical(sigma0_m = sigma_0_m, w = w)
#c_times = IHf.base_run.sh.t_collapse(sigma_eff_m = sigma_eff_array, rmax = rmax_array, Vmax = vmax_array) / sh.Gyr


    def universe_age_calc(self, zf_max = 10, zmin_N = 1000):
        """
        Computes the age of the universe using equation (2.22) from https://arxiv.org/pdf/2403.16633. However since we need this value to be computed at different
        redshifts we compute the age of the universe at different redshifts and return an interpolation function for this.
        Returns interpolation function for the age of the universe at a given redshift. The interpolation function returns the age of the universe in seconds.
        """
        self.zf_max = zf_max
        z_dummy_min = np.linspace(0, zf_max, zmin_N)
        z_dummy = np.linspace(z_dummy_min, 1000, 10000)
        t_U = simpson(1./(self.sh.Hubble(z_dummy)*(1.+z_dummy)),x=z_dummy,axis=0)
        self.universe_age = make_interp_spline(z_dummy_min, t_U)
        #return U_age

    def U_age(self, z):
        if not hasattr(self, "universe_age"):
            raise StageException("The universe_age_calc function needs to be called first before calling U_age")
        z_too_large = np.any(z > self.zf_max)
        if z_too_large:
            print(f"Warning, one (or more) redshift values given to the U_age function is too large and have been limited to {self.zf_max}")
        z = np.where(z > self.zf_max, self.zf_max, z)
        return self.universe_age(z)
        
        

    def compute_M0_from_Mz(self, masses, redshifts, use_fit = True):
        """
        Returns the masses of halos at z = 0 from the values of the mass at a redshift z, so (M(z)). The fit was done using simulated data where M_0 was computed from M(z) using the method described in https://arxiv.org/pdf/1409.5228. This same method is used by Sashimi
        """
        md = -0.122
        zd = 0.179
        mm = 1.019

        def M0_fit(z_mass, md, zd, mm):
            z, log10mass = z_mass
            return md + z * zd + mm * log10mass

        def M0_from_interpolation(self, redshifts, masses, N = 1000):
            """
            More accurate way to compute M0 based on M(z), although this will take much longer to compute compared to the fitted function.
            """
            M0_computed = []
            for z, mass in zip(redshifts, masses):
                M0_array = np.logspace(0, 3, N) * mass
                M0_interp = make_interp_spline(self.sh.Mzi(M0_array, z), M0_array)
                M0_computed.append(M0_interp(mass))
            M0_array = np.array(M0_computed)
            return M0_array

        def M0_from_fit(redshifts, masses):
            return M0_fit((redshifts, masses), md, zd, mm)

        if use_fit:
            return M0_from_fit(redshifts, masses)
        else:
            return M0_from_interpolation(self, redshifts, masses)

    def compute_z05_from_Mz(self, redshifts, masses, mfraction = 0.5, use_fit = True, use_M0_fit = True, precomputed_M0 = False, M0_array = None):
        zd = 2.63
        mz = 0.993
        mm = 0.137


        def simple_zfit(z_mass, zd, mz, mm):
            z, log10mass = z_mass
            return zd + mz*z - mm * log10mass


        def z_05_from_interpolation(self, M0_array, masses, N = 1000):
            z_05 = []
            z_dummy = np.linspace(0, 100, N)
            for M0, mass in zip(M0_array, masses):
                Mz_z = self.sh.Mzi(M0, z_dummy)
                Mz_interp = make_interp_spline(np.flip(Mz_z), np.flip(z_dummy))
                z_05.append(Mz_interp(mass * mfraction))
            z_05 = np.array(z_05)
            return z_05
            
        def z_05_from_fit(redshifts, masses):
            return simple_zfit((redshifts, masses), zd, mz, mm)

        if use_fit:
            return z_05_from_fit(redshifts, masses)
        elif precomputed_M0:
            return z_05_from_interpolation(self, M0_array, masses)
        else:
            return z_05_from_interpolation(self, self.compute_M0_from_Mz(redshifts, masses, use_M0_fit))
        

    @staticmethod
    def get_rmax(rs):
        return 2.163 * rs
        
    def get_vmax(self, rs, rhos):
        return rs * np.sqrt(4 * np.pi * self.sh.G * rhos / 4.625)


"""
seed_min, seed_max = 0, 1e9
first_seed = np.random.randint(low = seed_min, high = seed_max)
HC = IH.base_run

# Computing properties at z = 0
z0_array = np.zeros(len(M0_array))
c200_array_z0 = HC.concentration_mass_200(M0_array, z0_array)
r200_array_z0 = HC.get_r200_from_m200(M0_array, z0_array)
np.random.seed(first_seed)
c200_array_z0 = 10 ** np.random.normal(np.log10(c200_array_z0), HC.IHB.sigmalogc)
rs_array_z0 = r200_array_z0 / c200_array_z0
rhos_array_z0 = M0_array / (4 * np.pi * rs_array_z0 ** 3 * HC.f_NFW(c200_array_z0))

c200_array_z = HC.concentration_mass_200(masses, redshifts)
r200_array_z = HC.get_r200_from_m200(masses, redshifts)
np.random.seed(first_seed)
c200_array_z = 10 ** np.random.normal(np.log10(c200_array_z), HC.IHB.sigmalogc)
rs_array_z = r200_array_z / c200_array_z
rhos_array_z = masses / (4 * np.pi * rs_array_z ** 3 * HC.f_NFW(c200_array_z))

# Computing collapse times 
IHB = IH.base_run.IHB
sigma_eff_interp = sh.sigma_eff_m_interpolate_analytical(sigma0_m = IHB.sigma_0 * sh.cm **2 / sh.gram,
                                                                   w = IHB.vel_scale * sh.km / sh.s)
collapse_time = sh.t_collapse
def get_rmax(rs):
    return 2.163 * rs

def get_vmax(rs, rhos):
    return rs * np.sqrt(4 * np.pi * sh.G * rhos / 4.625)
rmax_0 = get_rmax(rs_array_z0 * sh.kpc)
vmax_0 = get_vmax(rs_array_z0 * sh.kpc, rhos_array_z0 * sh.Msun / (sh.kpc**3))
sigma_eff_array = sigma_eff_interp(vmax_0)
tc = sh.t_collapse(sigma_eff_m = sigma_eff_array, rmax = rmax_0, Vmax = vmax_0) / sh.Gyr
halo_ages = (U_age(redshifts) - U_age(field_halo_fz * np.ones(len(redshifts)))) / sh.Gyr
time_ratio = halo_ages / tc
#sh_tt = IH.base_run.SH.S_tt_ratio
sh_zacc = IH.base_run.SH.S_z_acc
vmax_sub = IH.base_run.SH.S_VmaxCDM_z0
rs_sub = IH.base_run.SH.S_rsCDM_z0
rhos_sub = IH.base_run.SH.S_rhosCDM_z0
rmax_sub = IH.base_run.SH.S_rmaxCDM_z0
sigma_eff_m = sigma_eff_interp(vmax_sub)
ctimes_sub = IH.base_run.sh.t_collapse(sigma_eff_m = sigma_eff_m, rmax = rmax_sub, Vmax = vmax_sub) / sh.Gyr



# Computing SIDM properties from properties at z = 0 and tc:
rs_array_z0
rhos_array_z0
time_ratio

get_rhos = sh.get_rhos
get_rs = sh.get_rs
get_rc = sh.get_rc

time_ratio_11 = np.where(time_ratio > 1.1, 1.1, time_ratio)
collapsed = np.where(time_ratio >= 1.0)


rhos_array_z = get_rhos(rhos_array_z0, time_ratio_11)
rs_array_z = get_rs(rs_array_z0, time_ratio_11)
rc_array_z = get_rc(rs_array_z0, time_ratio_11)
"""




















