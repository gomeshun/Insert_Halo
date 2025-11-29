from pyHalo.Halos.HaloModels.TNFW import TNFWSubhalo
from pyHalo.Halos.HaloModels.TNFWFromParams import TNFWFromParams
from IH_errors import StageException
import numpy as np


class TNFWSubhaloSASHIMI(TNFWSubhalo):
    def __init__(self, mass, x, y, r3d, z,
                 sub_flag, lens_cosmo_instance, args,
                 truncation_class, concentration_class, unique_tag):
        TNFWSubhalo.__init__(self, mass, x, y, r3d, z,
                 sub_flag, lens_cosmo_instance, args,
                 truncation_class, concentration_class, unique_tag)

class TNFWSubHaloSASHIMI_FP(TNFWFromParams):
    def __init__(self, mass, x_kpc, y_kpc, r3d, z,sub_flag,
                 lens_cosmo_instance, args, unique_tag = None):
        TNFWFromParams.__init__(self, mass, x_kpc, y_kpc, r3d, z,sub_flag,
                 lens_cosmo_instance, args, unique_tag)

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

                    zlens = self.SDG.IHB.zlens
                    arcsec_per_kpc = self.cosmology.LCDM_flat.arcsec_per_kpc_proper(zlens).value

                    lens_cosmo_instance = self.PH.cdm.lens_cosmo

                    if 1 - np.sum(survive)/len(survive) > 0.05:
                        print("Warning! The fraction of tidally disrupted subhalos compared to the total number of subhalos is larger than 5%")
                    

                    if len(m_z0) != len(x_arcsec):
                        raise IndexError("The selected Sashimi subhalo data should be of the same length as the generated position data in generate_new_halo_classes")
                    for mass, x, y, rt, rs, rhos, z_infall, sur in zip(m_z0, x_arcsec, y_arcsec, rt_z0, rs_z0, rhos_z0, z_acc, survive):
                        if sur == 1:
                            rand = np.random.rand()
                            args = {"r_trunc_kpc": rt*rs, "rs": rs, "rhos": rhos,
                                   "rv": self.cosmology.get_virial_radius(mass, zlens), "z_infall": z_infall,
                                   "index": rand}
                            self.new_halos.append(TNFWSubHaloSASHIMI_FP(mass, x / arcsec_per_kpc, y / arcsec_per_kpc, None, zlens, True, 
                                                                    lens_cosmo_instance, args, rand))
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


























