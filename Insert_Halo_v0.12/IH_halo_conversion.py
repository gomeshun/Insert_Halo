from pyHalo.Halos.HaloModels.TNFW import TNFWSubhalo
from pyHalo.Halos.HaloModels.TNFWFromParams import TNFWFromParams
from IH_errors import StageException
import numpy as np
import scipy.integrate as integrate


class TNFWSubhaloSASHIMI(TNFWSubhalo):
    def __init__(self, mass, x, y, r3d, z,
                 sub_flag, lens_cosmo_instance, args,
                 truncation_class, concentration_class, unique_tag):
        TNFWSubhalo.__init__(self, mass, x, y, r3d, z,
                 sub_flag, lens_cosmo_instance, args,
                 truncation_class, concentration_class, unique_tag)


class TNFWSubHaloSASHIMI_FP(TNFWFromParams):
    def __init__(self, mass, x_kpc, y_kpc, r3d, z,sub_flag,
                 lens_cosmo_instance, args, unique_tag = None, r200 = 0):
        TNFWFromParams.__init__(self, mass, x_kpc, y_kpc, r3d, z,sub_flag,
                 lens_cosmo_instance, args, unique_tag)
        self.r200 = r200

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

                    zlens = self.SDG.IHB.zlens
                    arcsec_per_kpc = self.cosmology.LCDM_flat.arcsec_per_kpc_proper(zlens).value

                    lens_cosmo_instance = self.PH.cdm.lens_cosmo

                    if 1 - np.sum(survive)/len(survive) > 0.05:
                        print("Warning! The fraction of tidally disrupted subhalos compared to the total number of subhalos is larger than 5%")
                    

                    if len(m_z0) != len(x_arcsec):
                        raise IndexError("The selected Sashimi subhalo data should be of the same length as the generated position data in generate_new_halo_classes")
                    for mass, x, y, rt, rs, rhos, z_infall, sur, m200_acc in zip(m_z0, x_arcsec, y_arcsec, rt_z0, rs_z0, rhos_z0, z_acc, survive, m200):
                        if sur == 1:
                            rand = np.random.rand()
                            args = {"r_trunc_kpc": rt*rs, "rs": rs, "rhos": rhos,
                                   #"rv": self.cosmology.concentration_calc(mass , zlens) * rs, "z_infall": z_infall,
                                   #"rv": self.cosmology.get_virial_radius(m200_acc, z_infall), "z_infall": z_infall,
                                   #"rv": rt * rs, "z_infall": z_infall,
                                   "rv": self.compute_Rvir_from_M200(m200_acc, z_infall, rs),
                                   "index": rand, "z_infall": z_infall}
                            self.new_halos.append(TNFWSubHaloSASHIMI_FP(mass, x / arcsec_per_kpc, y / arcsec_per_kpc, None, zlens, True, 
                                                                    lens_cosmo_instance, args, rand, r200 = self.cosmology.get_virial_radius(m200_acc, z_infall)))
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
        a1 = 0.5116
        a2 = -0.4283
        a3 = -3.13e-3
        a4 = -3.52e-5
        p = a2 + a3 * np.log(f) + a4 * (np.log(f)) ** 2
        return (a1 * f ** (2 * p) + (3 / 4) ** 2) ** (-1/2) + 2 * f

    def f_x(self, x):
        return x ** 3 * (np.log(1 + x ** -1) - (1 + x) ** -1)

    def compute_Rvir_from_M200(self, m200, z_acc, rs):
        c200 = self.get_r200_from_m200(m200, z_acc) / rs
        Delta_vir = self.Delta_vir(z_acc)
        rs_rvir = self.x_f(Delta_vir / 200 * self.f_x(1 / c200))
        #Mvir_M200 = Delta_vir / 200 / rs_rvir / c200
        #return Mvir_M200 * m200
        return rs / rs_rvir
        











