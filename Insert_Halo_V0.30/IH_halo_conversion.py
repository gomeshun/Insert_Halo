from pyHalo.Halos.HaloModels.TNFW import TNFWSubhalo
from pyHalo.Halos.HaloModels.TNFWFromParams import TNFWFromParams
from IH_errors import StageException
import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import make_interp_spline


class TNFWSubhaloSASHIMI(TNFWSubhalo):
    def __init__(self, mass, x, y, r3d, z,
                 sub_flag, lens_cosmo_instance, args,
                 truncation_class, concentration_class, unique_tag):
        TNFWSubhalo.__init__(self, mass, x, y, r3d, z,
                 sub_flag, lens_cosmo_instance, args,
                 truncation_class, concentration_class, unique_tag)


class TNFWSubHaloSASHIMI_FP(TNFWFromParams):
    def __init__(self, mass, x_kpc, y_kpc, r3d, z,sub_flag,
                 lens_cosmo_instance, args, unique_tag = None, r200 = 0, rt_sh = 0, mass_evolved = None):
        TNFWFromParams.__init__(self, mass, x_kpc, y_kpc, r3d, z,sub_flag,
                 lens_cosmo_instance, args, unique_tag)
        self.r200 = r200
        self.rt_sh = rt_sh
        self.mass_evolved = mass_evolved

    def dens_TNFW(self, r):
        args = self.params_physical
        rs = args["rs"]
        rt = args["r_trunc_kpc"]
        rhos = args["rhos"]
        x = r / rs
        tau = rt / rs
        return rhos * tau**2 / ((x * (x + 1)**2) * (tau**2 + x**2))




class HaloConversion():
    def __init__(self, spatial_distribution_generation):
        self.SDG = spatial_distribution_generation
        self.SH = self.SDG.SH
        self.PH = self.SDG.PH
        self.cosmology = self.SDG.cosmology
        self.from_params = True
        self.new_halos = []
        self.sh = self.SH.sh_settings
        self.GM = self.SDG.GM
        self.IHB = self.SDG.IHB



    def generate_new_halo_classes(self):
        if self.from_params:
            if not self.SH.ran_sashimi:
                raise StageException("generate_new_halo_classes requires that Sashimi was run first.")
            else:
                """
                randomization = self.SDG.IHB.randomization
                if not randomization and self.SDG.Nsub != None:
                    Nsub = self.SDG.Nsub
                elif randomization and self.SDG.Nsub_random != None:
                    Nsub = self.SDG.Nsub_random
                else:
                    if randomization:
                        Nsub = self.SDG.rand_Nsub()
                    else:
                        Nsub = self.SDG.calc_Nsub()
                """
                if self.SDG.IHB.DM_type == "CDM":
                    z_acc = self.SH.S_z_acc
                    m_z0 = self.SH.S_m_z0 / self.sh.Msun
                    rs_z0 = self.SH.S_rsCDM_z0 / self.sh.kpc
                    rhos_z0 = self.SH.S_rhosCDM_z0 / self.sh.Msun * self.sh.kpc**3
                    rt_z0 = self.SH.S_ctCDM_z0
                    survive = self.SH.S_surviveCDM
                    x_arcsec = self.SDG.x_arcsec
                    y_arcsec = self.SDG.y_arcsec

                    m200 = self.SH.S_ma200 / self.sh.Msun

                    rs_acc = self.SH.S_rsCDM_acc / self.sh.kpc
                    rhos_acc = self.SH.S_rhosCDM_acc / self.sh.Msun * self.sh.kpc**3
                    

                    zlens = self.SDG.IHB.zlens
                    arcsec_per_kpc = self.cosmology.LCDM_flat.arcsec_per_kpc_proper(zlens).value

                    lens_cosmo_instance = self.PH.cdm.lens_cosmo

                    if 1 - np.sum(survive)/len(survive) > 0.05:
                        print("Warning! The fraction of tidally disrupted subhalos compared to the total number of subhalos is larger than 5%")
                    

                    if len(m_z0) != len(x_arcsec) or len(m_z0) != len(y_arcsec):
                        raise IndexError("The selected Sashimi subhalo data should be of the same length as the generated position data in generate_new_halo_classes")
                    for mass, x, y, rt, rs, rs_a, rhos, z_infall, sur, m200_acc in zip(m_z0, x_arcsec, y_arcsec, rt_z0, rs_z0, rs_acc, rhos_z0, z_acc, survive, m200):
                        if sur == 1:
                            rand = np.random.rand()
                            args = {"r_trunc_kpc": self.find_new_rt(rt * rs, rs, rhos, mass), "rs": rs_a, "rhos": rhos,
                                   #{"r_trunc_kpc": rt*rs, "rs": rs, "rhos": rhos,
                                   
                                   
                                   #{"r_trunc_kpc": 2 * self.cosmology.get_virial_radius(m200_acc, z_infall), "rs": rs, "rhos": rhos,
                                   #"rv": self.cosmology.concentration_calc(mass , zlens) * rs, "z_infall": z_infall,
                                   #"rv": self.cosmology.get_virial_radius(m200_acc, z_infall), "z_infall": z_infall,
                                   #"rv": rt * rs, "z_infall": z_infall,
                                   "rv": self.compute_Rvir_from_M200(m200_acc, z_infall, rs_a),
                                   "index": rand, "z_infall": z_infall}

                            # First we initialize the TNFWSubHaloSASHIMI_FP class using the accretion rs values as the constructor of the 
                            # pyHalo class TNFWSubHaloFromParams computes the concentration using c = r_vir / rs, but these values are at 
                            # the moment of accretion. After initializing (and thus computing the concentration) we replace the rs value
                            # with the Sashimi rs value after tidal effects. This is possible since pyHalo assumes that rs remains constant
                            # over time, but rho_s does evolve.
                            new_subhalo = TNFWSubHaloSASHIMI_FP(m200_acc, x / arcsec_per_kpc, y / arcsec_per_kpc, None, zlens, True, 
                                                                    lens_cosmo_instance, args, rand, r200 = self.cosmology.get_virial_radius(m200_acc, z_infall),
                                                                       rt_sh = rt * rs, mass_evolved = mass)
                            new_subhalo.params_physical["rs"] = rs
                            self.new_halos.append(new_subhalo)
        return self.new_halos

    

    def insert_CDM_halos(self):
        if not self.PH.ran_pyhalo:
            raise StageException("pyhalo has to have been run in order to insert the new halos using the insert_CDM_halos function.")
        else:
            if len(self.PH.cdm.subhalos) > 0:
                index_new = 0
                for index, halo in enumerate(self.PH.cdm.halos):
                    if halo.is_subhalo:
                        self.PH.cdm.halos[index] = self.new_halos[index_new]
                        index_new += 1
                for halo in self.new_halos[index_new:]:
                    self.PH.cdm.halos.append(halo)
            else:
                self.PH.cdm.halos = self.PH.cdm.halos + self.new_halos

        self.PH.cdm._reset() #Resets all class attributes to the current set of halos contained within the realization
        # This reset function has to be run in order for the pyHalo realization to recognize the new halos inserted into it.

    def insert_CDM_halos_corr(self):
        """
        This function inserts the halos generated from Sashimi data into the pyHalo realization. This is the updated function to 
        insert_CDM_halos and this function should be used instead of that one.
        """
        # Checking whether pyHalo has run, without this we cannot insert halo into a pyHalo realization:
        if not self.PH.ran_pyhalo:
            raise StageException("pyhalo has to have been run in order to insert the new halos using the insert_CDM_halos function.")
        else:
            # Checking whether pyHalo generated some subhalos of its own:
            if len(self.PH.cdm.subhalos) > 0:

                # Removing pyHalo subhalos and replacing them with Sashimi subhalos:
                ph_indices = []
                for index, halo in enumerate(self.PH.cdm.halos):
                    if not halo.is_subhalo:
                        ph_indices.append(index)
                self.PH.cdm.halos = np.array(self.PH.cdm.halos)
                self.PH.cdm.halos = self.PH.cdm.halos[ph_indices]
                self.PH.cdm.halos = list(self.PH.cdm.halos)
                self.PH.cdm.halos += self.new_halos
            else:
                self.PH.cdm.halos += self.new_halos

        self.PH.cdm._reset() #Resets all class attributes to the current set of halos contained within the realization
        # This reset function has to be run in order for the pyHalo realization to recognize the new halos inserted into it.


    


    def avg_rho_all_subhalos(self):
        subhalos = self.PH.cdm.subhalos

        avg_densities = []

        for halo in subhalos:
            parameters = halo.params_physical
            avg_densities.append(self.cosmology.calc_avg_density(self.cosmology.SH_TNFW, R_max = parameters["rv"], 
                                                                 rs = parameters["rs"], rhos = parameters["rhos"],
                                                                 rt = parameters["r_trunc_kpc"]))

        return avg_densities


    #def get_r200_CDM(self, density_profile, rs, rhos, rt, R_max):
        
        
            
        
    def get_r200_from_m200(self, m200, z_acc):
        Delta_200 = 200
        rho_crit_z = self.sh.rhocrit(z_acc) / self.sh.Msun * self.sh.kpc**3
        r_200 = (3 * m200 / (4 * np.pi * Delta_200 * rho_crit_z)) ** (1/3)
        return r_200

    @staticmethod
    def c200(r200, rs):
        return r200 / rs


    def Delta_vir(self, z):
        x = self.sh.OmegaM * (1 + z) ** 3 / (self.sh.OmegaM * (1 + z) ** 3 + 1 - self.sh.OmegaM)
        return 18 * np.pi ** 2 + 82 * x - 39 * x ** 2


    def x_f(self, f):
        """
        Fitted function representing the inverted f(x) function, which is the dimensionless integrated NFW density
        https://arxiv.org/pdf/astro-ph/0203169
        """
        a1 = 0.5116
        a2 = -0.4283
        a3 = -3.13e-3
        a4 = -3.52e-5
        p = a2 + a3 * np.log(f) + a4 * (np.log(f)) ** 2
        return (a1 * f ** (2 * p) + (3 / 4) ** 2) ** (-1/2) + 2 * f

    def f_x(self, x):
        """
        Dimensionless integrated NFW density
        """
        return x ** 3 * (np.log(1 + x ** -1) - (1 + x) ** -1)

    def compute_Rvir_from_M200(self, m200, z_acc, rs):
        """
        Procedure found in the Appendix of https://arxiv.org/pdf/astro-ph/0203169 to convert M200 into Mvir, here we 
        stop a step before to get the virial radius r_vir instead.
        """
        c200 = self.get_r200_from_m200(m200, z_acc) / rs
        Delta_vir = self.Delta_vir(z_acc)
        rs_rvir = self.x_f(Delta_vir / 200 * self.f_x(1 / c200))
        #Mvir_M200 = Delta_vir / 200 / rs_rvir / c200
        #return Mvir_M200 * m200
        return rs / rs_rvir
        

    @staticmethod
    def TNFW_int_inf(rt, rs):
        """
        Expression for the total mass of the smoothly truncated NFW profile (also known as BMO profile or TNFW profile).
        The value returned is dimensionless as this expression uses tau = r_t / r_s. This equation is taken from
        https://arxiv.org/pdf/0705.0682
        """
        tau = np.where(rt > 0, rt / rs, np.max(rt/rs) / len(rt) / 10)
        frac = tau ** 2 / (tau ** 2 + 1) ** 2
        tail = np.pi * tau - tau ** 2 - 1 + (tau ** 2 - 1) * np.log(tau)
        return frac * tail

    @staticmethod
    def find_new_rt(rt_sh, rs, rhos, mass, N = 1000):
        """
        Estimates a new value for the truncation radius based on matching the total mass between the sharply truncated NFW profile and the 
        smoothly truncated NFW profile (or BMO or TNFW profile). This function uses the truncation radius from Sashimi (after multiplying the 
        Sashimi output by r_s to make it a physical distance)
        """
        rt_array = np.linspace(rt_sh / N / 10, rt_sh * 5, N)
        mass_comp = mass / (4 * np.pi * rhos * rs ** 3)
        dl_mass_int = HaloConversion.TNFW_int_inf(rt_array, rs)
        rt_spline_inv = make_interp_spline(dl_mass_int, rt_array)
        return rt_spline_inv(mass_comp)
        





