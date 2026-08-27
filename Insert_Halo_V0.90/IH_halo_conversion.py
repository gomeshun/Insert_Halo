#from pyHalo.Halos.HaloModels.TNFW import TNFWSubhalo
from pyHalo.Halos.HaloModels.TNFWFromParams import TNFWFromParams
from pyHalo.Halos.HaloModels.TNFW import TNFWFieldHalo
from IH_errors import StageException
from IH_errors import SettingsException
import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import make_interp_spline
from scipy.integrate import simpson
from scipy.integrate import trapezoid
from IH_alternative_DM_dist import SIDM_halo


class TNFWSubHaloSASHIMI_FP(TNFWFromParams):
    def __init__(self, mass, x_kpc, y_kpc, r3d, z,sub_flag,
                 lens_cosmo_instance, args, unique_tag = None, r200 = 0, rt_sh = 0, mass_evolved = None, z_acc = None):
        TNFWFromParams.__init__(self, mass, x_kpc, y_kpc, r3d, z,sub_flag,
                 lens_cosmo_instance, args, unique_tag)
        self.r200 = r200
        self.rt_sh = rt_sh
        self.mass_evolved = mass_evolved
        self.z_acc = z_acc

    def dens_TNFW(self, r):
        args = self.params_physical
        rs = args["rs"]
        rt = args["r_trunc_kpc"]
        rhos = args["rhos"]
        x = r / rs
        tau = rt / rs
        return rhos * tau**2 / ((x * (x + 1)**2) * (tau**2 + x**2))

    @property
    def rs(self):
        params = self.params_physical
        return params["rs"]

    @property
    def rhos(self):
        params = self.params_physical
        return params["rhos"]

    @property
    def rt(self):
        params = self.params_physical
        return params["r_trunc_kpc"]

    @property
    def profile_args(self):
        return 


class TNFWFieldHaloSASHIMI(TNFWFieldHalo):
    def __init__(self, mass, x, y, r3d, z,
                 sub_flag, lens_cosmo_instance, args,
                 truncation_class, concentration_class, unique_tag, r200):
        TNFWFieldHalo.__init__(self, mass, x, y, r3d, z,
                 sub_flag, lens_cosmo_instance, args,
                 truncation_class, concentration_class, unique_tag)
        self.r200 = r200




class HaloConversion():
    def __init__(self, spatial_distribution_generation):
        self.SDG = spatial_distribution_generation
        self.SH = self.SDG.SH
        self.PH = self.SDG.PH
        self.cosmology = self.SDG.cosmology
        self.new_subhalos = []
        self.new_fieldhalos = []
        self.sh = self.SH.sh_settings
        self.GM = self.SDG.GM
        self.IHB = self.SDG.IHB
        self.rt_subhalos_CDM = []
        self.rt_fieldhalos_CDM = []
        self.rhos_fieldhalos_NFW = []



    def generate_new_halo_classes(self):
        #======================================================
        #  CDM Subhalo class generation
        #======================================================
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

                lens_cosmo_instance = self.PH.DM_dist.lens_cosmo

                if 1 - np.sum(survive)/len(survive) > 0.05:
                    print("Warning! The fraction of tidally disrupted subhalos compared to the total number of subhalos is larger than 5%")
                

                if len(m_z0) != len(x_arcsec) or len(m_z0) != len(y_arcsec):
                    raise IndexError("The selected Sashimi subhalo data should be of the same length as the generated position data in generate_new_halo_classes")
                for mass, x, y, rt, rs, rs_a, rhos, z_infall, sur, m200_acc in zip(m_z0, x_arcsec, y_arcsec, rt_z0, rs_z0, rs_acc, rhos_z0, z_acc, survive, m200):
                    if sur == 1:
                        rand = np.random.rand()
                        #args = {"r_trunc_kpc": self.find_new_rt(rt * rs, rs, rhos, mass), "rs": rs_a, "rhos": rhos,
                               #{"r_trunc_kpc": rt*rs, "rs": rs, "rhos": rhos,
                        args = {"r_trunc_kpc": self.new_trunc_analytical(rs, rhos, mass), "rs": rs_a, "rhos": rhos,
                        #profile_properties = {"profile": "TNFW", "rs": rs, "rhos": rhos}
                        #rt_new = self.num_new_truncation(mass, 3, profile_properties, Nr = 100)
                        #self.rt_subhalos_CDM.append(rt_new)
                        #args = {"r_trunc_kpc": rt_new, "rs": rs_a, "rhos": rhos,
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
                                                                   rt_sh = rt * rs, mass_evolved = mass, z_acc = z_infall)
                        new_subhalo.params_physical["rs"] = rs
                        self.new_subhalos.append(new_subhalo)

                #======================================================
                #  CDM Fieldhalo class generation
                #======================================================
                if not self.PH.ran_pyhalo:
                    raise StageException("pyHalo has to have been run in order to generate new field halos from the pyHalo field halos in generate_new_halo_classes.")
                else:
                    field_halos = self.PH.DM_dist.field_halos
                    if len(field_halos) <= 0:
                        raise SettingsException("Current configuration contains no field halos so generate_new_halo_classes cannot be run fully.")
                    fhalo_masses = []
                    fhalo_z = []
                    fhalo_x = []
                    fhalo_y = []
        
                    for fhalo in field_halos:
                        fhalo_masses.append(fhalo.mass)
                        fhalo_z.append(fhalo.z)
                        fhalo_x.append(fhalo.x)
                        fhalo_y.append(fhalo.y)
                    
                    fhalo_masses = np.array(fhalo_masses)
                    fhalo_z = np.array(fhalo_z)
                    
                    c200_array = self.concentration_mass_200(fhalo_masses, fhalo_z)
                    r200_array = self.get_r200_from_m200(fhalo_masses, fhalo_z)

                    rt_array = r200_array * 2.0
                    
                    #c_scatter = 0.128
                    c200_array = 10 ** np.random.normal(np.log10(c200_array), self.IHB.sigmalogc)
                    
                    rs_array = r200_array / c200_array
                    rhos_array_NFW = fhalo_masses / (4 * np.pi * rs_array ** 3 * self.f_NFW(c200_array))
                    rhos_array = fhalo_masses / (4 * np.pi * rs_array ** 3 * self.ftau_x(rt_array / rs_array, r200_array / rs_array))
                
                    #print(np.min(rhos_array), np.max(rhos_array))
        
                    lens_cosmo_instance = self.PH.DM_dist.lens_cosmo
        
                    # Generating the pyHalo field halos based on generated masses and redshifts:
                    for mass, z, x, y, c200, rs, rhos, r200, rhos_NFW in zip(fhalo_masses, fhalo_z, fhalo_x, fhalo_y, c200_array, rs_array, rhos_array, r200_array, rhos_array_NFW):
        
                        args = {}
                        rand = np.random.rand()
                        new_fieldhalo = TNFWFieldHaloSASHIMI(mass, x, y, None, z, False, lens_cosmo_instance, args, None, None, rand, r200)
                        new_fieldhalo._c = self.compute_Rvir_from_M200(mass, z, rs) / rs #Converting first to get concentration for the mass inclosed in r_vir instead of r_200
                        #rt = self.find_new_rt(rt_sh, rs, rhos, mass, N = 1000)
        
                        # Finding the fieldhalo truncation radius:
                        tau_temp = np.linspace(0.001, 1000, 1000) # Careful to make sure this does not start at '0'.
                        tau_inv = make_interp_spline(self.ftau_x(tau_temp, c200), tau_temp)
                        tau_new = tau_inv(mass / (4 * np.pi * rs**3 * rhos)).tolist()
                        rt_new = tau_new * rs
        
                        #tau_inv = make_interp_spline(self.ftau_inf(tau_temp), tau_temp)
                        #tau_new = tau_inv(mass / (4 * np.pi * rs **3 * rhos))
                        #profile_properties = {"profile": "TNFW", "rs": rs, "rhos": rhos}
                        #print(profile_properties["profile"])
                        #if mass > 1e9:
                        #    rt_new = self.num_new_truncation(mass, np.log10(r200), profile_properties, Nr = 1000)
                        #else:
                        #    rt_new = self.num_new_truncation(mass, np.log10(r200), profile_properties)
                        #rt_new = r200 * 2
                        new_fieldhalo._profile_args = (new_fieldhalo._c, rt_new)
                        new_fieldhalo._params_physical = {"rhos": rhos, "rs": rs, "r200": r200, "r_trunc_kpc": rt_new}
                        new_fieldhalo._nfw_params = [rhos, rs, r200]
                        self.rt_fieldhalos_CDM.append(rt_new)
                        self.new_fieldhalos.append(new_fieldhalo)
                        self.rhos_fieldhalos_NFW.append(rhos_NFW)

            #======================================================
            #   SIDM subhalo generation
            #======================================================
            elif self.SDG.IHB.DM_type == "SIDM":
                z_acc = self.SH.S_z_acc
                m_z0 = self.SH.S_m_z0 / self.sh.Msun
                rs_z0 = self.SH.S_rsSIDM_z0 / self.sh.kpc
                rhos_z0 = self.SH.S_rhosSIDM_z0 / self.sh.Msun * self.sh.kpc**3
                rt_z0 = self.SH.S_ctCDM_z0
                survive = self.SH.S_surviveSIDM
                rc_z0 = self.SH.S_rcSIDM_z0 / self.sh.kpc
                x_arcsec = self.SDG.x_arcsec
                y_arcsec = self.SDG.y_arcsec
                t_over_tc = self.SH.S_tt_ratio

                m200 = self.SH.S_ma200 / self.sh.Msun

                rs_acc = self.SH.S_rsCDM_acc / self.sh.kpc
                rhos_acc = self.SH.S_rhosCDM_acc / self.sh.Msun * self.sh.kpc**3
                
                #print(np.max(t_over_tc))
                zlens = self.SDG.IHB.zlens
                #arcsec_per_kpc = self.cosmology.LCDM_flat.arcsec_per_kpc_proper(zlens).value
                #kpc_per_arcsec = self.GM.kpc_per_arcsec()

                #lens_cosmo_instance = self.PH.DM_dist.lens_cosmo

                if 1 - np.sum(survive)/len(survive) > 0.05:
                    print("Warning! The fraction of tidally disrupted subhalos compared to the total number of subhalos is larger than 5%")

                if len(m_z0) != len(x_arcsec) or len(m_z0) != len(y_arcsec):
                    raise IndexError("The selected Sashimi subhalo data should be of the same length as the generated position data in generate_new_halo_classes")
                for mass, x, y, rt, rs, rs_a, rhos, z_infall, sur, rc, t_tc in zip(m_z0, x_arcsec, y_arcsec, rt_z0, rs_z0, rs_acc, rhos_z0, z_acc, survive, rc_z0, t_over_tc):
                    if sur == 1:
                        profile_properties = {"profile": "TNFWC_SH_tau", "rhos": rhos, "rs": rs, "rc": rc}
                        args = profile_properties
                        #args["rt"] = self.num_new_truncation(mass, 3, properties)
                          
                        if t_tc > 1:
                            has_collapsed = True
                        else:
                            has_collapsed = False

                        if has_collapsed:
                            args["rt"] = self.num_new_truncation(mass, 3, profile_properties, Nr = 1000, log_rmin = -3, rt_logmin = np.log10(rc) - 3, Nt = 1000)

                        else:
                            rt_logmin, rt_logmax = -3, 2
                            args["rt"] = self.num_new_truncation(mass, 3, profile_properties, Nr = 1000, rt_logmin = rt_logmin, rt_logmax = rt_logmax)

                        if args["rt"] > 10**rt_logmax:
                            new_subhalo.mass_status = "incorrect"
                        
                        new_subhalo = SIDM_halo(mass, x, y, zlens, True, args, has_collapsed, t_tc)
                        #self.PH.SIDM_dist.halo_list.append(new_subhalo)
                        new_subhalo.z_acc = z_infall
                        self.new_subhalos.append(new_subhalo)


                #======================================================
                #   SIDM fieldhalo generation
                #======================================================
                if not self.PH.ran_pyhalo:
                    raise StageException("pyhalo has to have been run in order to generate new field halos from the pyHalo field halos in generate_new_halo_classes.")
                else:
                    field_halos = self.PH.DM_dist.field_halos
                    if len(field_halos) <= 0:
                        raise SettingsException("Current configuration contains no field halos so generate_new_halo_classes cannot be run fully.")
                    fhalo_masses = []
                    fhalo_z = []
                    fhalo_x = []
                    fhalo_y = []
        
                    for fhalo in field_halos:
                        fhalo_masses.append(fhalo.mass)
                        fhalo_z.append(fhalo.z)
                        fhalo_x.append(fhalo.x)
                        fhalo_y.append(fhalo.y)
                    
                    fhalo_masses = np.array(fhalo_masses)
                    fhalo_z = np.array(fhalo_z)
                    fhalo_x = np.array(fhalo_x)
                    fhalo_y = np.array(fhalo_y)

                    SIDM_calc = self.PH.SIDM_calc

                    # Generate random seed (for concentration)
                    seed_min, seed_max = 0, 1e9
                    seed_ = np.random.randint(low = seed_min, high = seed_max)
                    
                    # Compute CDM rs and rhos at z 
                    c200_array_z = self.concentration_mass_200(fhalo_masses, fhalo_z)
                    r200_array_z = self.get_r200_from_m200(fhalo_masses, fhalo_z)
                    np.random.seed(seed_)
                    c200_array_z = 10 ** np.random.normal(np.log10(c200_array_z), self.IHB.sigmalogc)
                    rs_array_z = r200_array_z / c200_array_z
                    rhos_array_z = fhalo_masses / (4 * np.pi * rs_array_z ** 3 * self.f_NFW(c200_array_z))
                    
                    # Compute M0_array
                    M0_array = SIDM_calc.compute_M0_from_Mz(fhalo_masses, fhalo_z, use_fit = False)

                    # Reset seed
                    np.random.seed(seed_)

                    # Compute rs and rhos at z = 0
                    z0_array = np.zeros(len(M0_array))
                    c200_array_z0 = self.concentration_mass_200(M0_array, z0_array)
                    r200_array_z0 = self.get_r200_from_m200(M0_array, z0_array)
                    c200_array_z0 = 10 ** np.random.normal(np.log10(c200_array_z0), self.IHB.sigmalogc)
                    rs_array_z0 = r200_array_z0 / c200_array_z0
                    rhos_array_z0 = M0_array / (4 * np.pi * rs_array_z0 ** 3 * self.f_NFW(c200_array_z0))

                    # Computing collapse times 
                    sigma_eff_interp = SIDM_calc.sigma_eff_interp
                    collapse_time = SIDM_calc.collapse_time
                    rmax_z = SIDM_calc.get_rmax(rs_array_z * self.sh.kpc)
                    vmax_z = SIDM_calc.get_vmax(rs_array_z * self.sh.kpc, rhos_array_z * self.sh.Msun / (self.sh.kpc**3))
                    sigma_eff_array = sigma_eff_interp(vmax_z)
                    tc = collapse_time(sigma_eff_m = sigma_eff_array, rmax = rmax_z, Vmax = vmax_z) / self.sh.Gyr

                    # Compute halo ages (and formation times)
                    SIDM_calc.universe_age_calc()
                    z05_array = SIDM_calc.compute_z05_from_Mz(fhalo_z, fhalo_masses, mfraction = self.IHB.formation_mass_fraction, use_fit = False, use_M0_fit = False, precomputed_M0 = True, M0_array = M0_array)
                    halo_ages = (SIDM_calc.U_age(fhalo_z) - SIDM_calc.U_age(z05_array)) / self.sh.Gyr
                    time_ratio = halo_ages / tc

                    # Compute SIDM properties using the basic method
                    time_ratio_lim1 = np.where(time_ratio > self.IHB.sh_tt_th, self.IHB.sh_tt_th, time_ratio)
                    
                    SIDMrhos_z = self.sh.get_rhos(rhos_array_z0, time_ratio_lim1)
                    SIDMrs_z = self.sh.get_rs(rs_array_z0, time_ratio_lim1)
                    SIDMrc_z = self.sh.get_rc(rs_array_z0, time_ratio_lim1)

                    #print(rs_array_z0)
                    #print(time_ratio_lim1)

                    #print(np.max(time_ratio))
                    # Generate SIDM_halo classes
                    for mass, x, y, rs, rhos, rc, t_tc, z, c200, rhosCDM, rsCDM, z05, t_c in zip(fhalo_masses, fhalo_x, fhalo_y, SIDMrs_z, SIDMrhos_z, SIDMrc_z, time_ratio, fhalo_z, c200_array_z, rhos_array_z, rs_array_z, z05_array, tc):
                        profile_properties = {"profile": "TNFWC_SH", "rhos": rhos, "rs": rs, "rc": rc}
                        args = profile_properties
                        #args["rt"] = self.num_new_truncation(mass, 3, properties)
                        #print(f"c200: {c200}")
                        #print(f"rs: {rs}")


                        #CDM_profile_properties = {"profile": "TNFW", "rhos": rhosCDM, "rs": rsCDM}

                        # Computing r_t for an equivalent CDM halo:
                        #if mass > 1e9:
                        #    CDM_profile_properties["rt"] = self.num_new_truncation(mass, np.log10(c200 * rsCDM), CDM_profile_properties, Nr = 1000)
                        #else:
                        #    CDM_profile_properties["rt"] = self.num_new_truncation(mass, np.log10(c200 * rsCDM), CDM_profile_properties)

                        # Then we compute the total mass of this CDM halo, which should be the same as the SIDM halo mass:
                        #r_array = np.logspace(-3, 4, int(np.log10(mass) * 1000))
                        #mass_total_CDM = self.mass_int_3d(r_array, CDM_profile_properties)

                        #BMO_dl_x = mass / (4*np.pi * rhos * rs**3)
                        #tau_temp = np.logspace(-2, 3, 1000)
                        #BMO_dl_match = self.ftau_x(tau_temp, c200)
                        #mass_to_rt = make_interp_spline(BMO_dl_match, tau_temp)
                        #au_CDM = mass_to_rt(BMO_dl_x)

                        #mass_total_CDMa = self.ftau_inf(tau_CDM) * 4*np.pi *rhos*rs**3
                        if t_tc >= 1.0:
                            has_collapsed = True
                        else:
                            has_collapsed = False                        

                        # Then compute the truncation radius for the SIDM profile, assuming the total mass of the SIDM halo is the same:

                        
                        #if mass > 1e9:
                        #    if has_collapsed: 
                        #       args["rt"] = self.num_new_truncation(mass, c200 * rsCDM, profile_properties, Nr = 10000, log_rmin = -4)
                        #    else:
                        #        args["rt"] = self.num_new_truncation(mass, c200 * rsCDM, profile_properties, Nr = 1000)
                        #lse:
                        #    if has_collapsed:
                        #        args["rt"] = self.num_new_truncation(mass, c200 * rsCDM, profile_properties, Nr = 1000)
                        #   else:
                        #        args["rt"] = self.num_new_truncation(mass, c200 * rsCDM, profile_properties, Nr = 100)
                                
                        #mass, x, y, z, sub_flag, profile_args, has_collapsed, t_over_tc
                        new_fieldhalo = SIDM_halo(mass, x, y, z, False, args, has_collapsed, t_tc)
                        #self.PH.SIDM_dist.halo_list.append(new_fieldhalo)
                        new_fieldhalo.r200 = c200 * rs
                        #new_fieldhalo.mass_total = {"masstn": mass_total_CDM, "massta": mass_total_CDMa}
                        #new_fieldhalo.rt_CDM = {"rta": tau_CDM * rs, "rtn": CDM_profile_properties["rt"]}
                        new_fieldhalo.z05 = z05
                        new_fieldhalo.tc = t_c
                        self.new_fieldhalos.append(new_fieldhalo)

                    # Add CDM component to each corresponding halo class
                    #if self.SIDM_keep_CDM:
                    #    for fhalo in self.new_fieldhalos:
                    #        fhalo.CDMproperties = {}
                        
                    
        return #self.new_halos

    

    def insert_CDM_halos(self):
        if not self.PH.ran_pyhalo:
            raise StageException("pyhalo has to have been run in order to insert the new halos using the insert_CDM_halos function.")
        else:
            if len(self.PH.DM_dist.subhalos) > 0:
                index_new = 0
                for index, halo in enumerate(self.PH.DM_dist.halos):
                    if halo.is_subhalo:
                        self.PH.DM_dist.halos[index] = self.new_halos[index_new]
                        index_new += 1
                for halo in self.new_halos[index_new:]:
                    self.PH.DM_dist.halos.append(halo)
            else:
                self.PH.DM_dist.halos = self.PH.DM_dist.halos + self.new_halos

        self.PH.DM_dist._reset() #Resets all class attributes to the current set of halos contained within the realization
        # This reset function has to be run in order for the pyHalo realization to properly recognize the new halos inserted into it.

    def insert_CDM_halos_corr(self):
        """
        This function inserts the halos generated from Sashimi data into the pyHalo realization. This is the updated function to 
        insert_CDM_halos and this function should be used instead of that one.
        """
        # Checking whether pyHalo has run, without this we cannot insert halo into a pyHalo realization:
        if not self.PH.ran_pyhalo:
            raise StageException("pyhalo has to have been run in order to insert the new halos using the insert_CDM_halos function.")
        else:
            assert len(self.PH.DM_dist.field_halos) == len(self.new_fieldhalos), "Number of new and existing field halos does not match in insert_CDM_halos_corr"
            # Checking whether pyHalo generated some subhalos of its own:
            if len(self.PH.DM_dist.subhalos) > 0:
                # Removing pyHalo subhalos and replacing them with Sashimi subhalos:
                ph_indices = []
                nfh_index = 0
                for index, halo in enumerate(self.PH.DM_dist.halos):
                    if not halo.is_subhalo:
                        ph_indices.append(index)
                        self.PH.DM_dist.halos[index] = self.new_fieldhalos[nfh_index]
                        nfh_index += 1
                self.PH.DM_dist.halos = np.array(self.PH.DM_dist.halos)
                self.PH.DM_dist.halos = self.PH.DM_dist.halos[ph_indices]
                self.PH.DM_dist.halos = list(self.PH.DM_dist.halos)
                self.PH.DM_dist.halos += self.new_subhalos
            else:
                self.PH.DM_dist.halos += self.new_subhalos
                nfh_index = 0
                for index, halo in enumerate(self.PH.DM_dist.halos):
                    if not halo.is_subhalo:
                        self.PH.DM_dist.halos[index] = self.new_fieldhalos[nfh_index]
                        nfh_index += 1
                

        self.PH.DM_dist._reset() #Resets all class attributes to the current set of halos contained within the realization
        # This reset function has to be run in order for the pyHalo realization to recognize the new halos inserted into it.

    def insert_SIDM_halos(self):
        # Inserting SIDM subhalos into alternative DM distribution class:
        for new_subhalo in self.new_subhalos:
            self.PH.SIDM_dist.halo_list.append(new_subhalo)
        for new_fieldhalo in self.new_fieldhalos:
            self.PH.SIDM_dist.halo_list.append(new_fieldhalo)
    


    def avg_rho_all_subhalos(self):
        subhalos = self.PH.DM_dist.subhalos

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

    def Omega(self, z):
        return (self.sh.OmegaM * (1 + z)**3) / (self.sh.OmegaM * (1 + z)**3 + self.sh.OmegaL) - 1


    def Delta_vir(self, z):
        #x = self.sh.OmegaM * (1 + z) ** 3 / (self.sh.OmegaM * (1 + z) ** 3 + 1 - self.sh.OmegaM)
        x = self.Omega(z)
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

    def compute_Rvir_from_M200(self, m200, z_acc, return_mvir = False):
        """
        Procedure found in the Appendix of https://arxiv.org/pdf/astro-ph/0203169 to convert M200 into Mvir, here we 
        stop a step before to get the virial radius r_vir instead.
        """
        c200 = self.concentration_mass_200(m200, z_acc)
        Delta_vir = self.Delta_vir(z_acc)
        rs_rvir = self.x_f(Delta_vir / 200 * self.f_x(1 / c200))
        r200 = self.get_r200_from_m200(m200, z_acc)
        #Mvir_M200 = Delta_vir / 200 / rs_rvir / c200
        #return Mvir_M200 * m200
        if return_mvir:
            Mvir_M200 = Delta_vir / 200 * (rs_rvir * c200)**(-3)
            return Mvir_M200 * m200
        else:
            rs = r200 / c200
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

    def new_trunc_analytical(self, rs, rhos, mass, N = 1000):
        tau_temp = np.logspace(-3, 3, N)
        ftau = self.ftau_inf(tau_temp)
        tau_ftau_interp = make_interp_spline(ftau, tau_temp)
        tau_found = tau_ftau_interp(mass / (4 * np.pi*rhos*rs**3))
        return tau_found * rs



    def concentration_mass_200(self, mass, z):
        """
        Concentration fitting function from https://arxiv.org/pdf/1502.00391, assuming the mass is given as m200.
        """
        logmass = np.log10(mass)
        alphazl4 = 1.7543 - 0.2766 * (1 + z) + 0.02039 * (1 + z) ** 2
        betazl4 = 0.2753 + 0.00351 * (1 + z) - 0.3038 * (1 + z) ** 0.0269
        gammazl4 = -0.01537 + 0.02102 * (1 + z) ** -0.1475
        c200_zl4 = 10 ** (alphazl4 + betazl4 * logmass * (1 + gammazl4 * (logmass) ** 2))
    
        alphazg4 = 1.3081 - 0.1078 * (1 + z) + 0.00398 * (1 + z) ** 2
        betazg4 = 0.0223 - 0.0944 * (1 + z) ** -0.3907
        c200_zg4 = 10 ** (alphazg4 + betazg4 * logmass)
    
        c200 = np.where(z > 4, c200_zg4, c200_zl4)
        return c200
    
    #def get_r200_from_m200(m200, z):
    #    Delta_200 = 200
    #    rho_crit_z = IHf.base_run.sh.rhocrit(z) / IHf.base_run.sh.Msun * IHf.base_run.sh.kpc**3
    #    r_200 = (3 * m200 / (4 * np.pi * Delta_200 * rho_crit_z)) ** (1/3)
    #    return r_200
    
    
    def f_NFW(self, c):
        return np.log1p(c) - c / (1 + c)


    def ftau_x(self, tau, x):
        """
        Full dimensionless integrated BMO profile (Appendix of https://arxiv.org/pdf/1101.0650)
        """
        first_term = tau**2 / (2 * (tau**2 + 1)**2 * (1 + x))
        second_term = -2 * (tau**2 + 1) * x + 4*tau * (x + 1) * np.arctan(x / tau) + (tau**2 - 1) * (1 + x) * np.log(tau**2 * (1 + x)**2 / (tau **2 + x **2))
        return first_term * second_term

    def ftau_inf(self, tau):
        """
        Same as ftau_x, except here the profile has been integrated until infinity (so it effectively equivalent to the total mass.
        """
        return tau**2 / (tau**2 + 1)**2 * (np.pi*tau - tau**2 - 1 + (tau**2 - 1)*np.log(tau))
  

        
    def mass_int_3d(self, r_array, properties):
        profile = getattr(self, properties["profile"])
        integrand = profile(r_array, properties) * r_array**2
        if r_array.ndim == 1:
            integral = 4 * np.pi * simpson(y = integrand * r_array, x = np.log(r_array))
        if r_array.ndim > 1:
            integral = 4 * np.pi * simpson(y = integrand * r_array, x = np.log(r_array), axis = 0) 
        
        #if r_array.ndim == 1:
        #    integral = 4 * np.pi * trapezoid(y = integrand * r_array, x = np.log(r_array))
        #elif r_array.ndim > 1:
        #    integral = 4 * np.pi * trapezoid(y = integrand * r_array, x = np.log(r_array), axis = 0)            
        return integral
    
    def num_new_truncation(self, mass, log_rmax, properties, Nr = 100, Nt = 100, log_rmin = -3, rt_logmin = -3, rt_logmax = 2):
        r_array = np.logspace(log_rmin, log_rmax, Nr)
        r_arrayd = np.expand_dims(r_array, axis = -1)
        r_array_2d = np.repeat(r_arrayd, Nt, axis = 1)
        rt_temp = np.logspace(rt_logmin, rt_logmax, Nt)
        properties["rt"] = rt_temp
        mass_array = self.mass_int_3d(r_array_2d, properties)
        #mass_nan = np.isnan(mass_array)
        #mass_notnan = np.invert(mass_nan)
        mass_inverted = make_interp_spline(mass_array, rt_temp)
        rt = mass_inverted(mass).tolist()
        return rt

    def TNFWC_SH_tau(self, r, properties):
        """
        Smoothly truncated TNFWC profile as used by https://arxiv.org/pdf/2403.16633 (expect this paper does not use the truncation term)
        """
        tau = properties["rt"] / properties["rs"]
        x = r / properties["rs"]
        beta = properties["rc"] / properties["rs"]
        rho = properties["rhos"] / ((x ** 4 + beta ** 4) ** (1/4) * (x + 1) ** 2) * tau**2 / (x **2 + tau**2)
        return rho

    def TNFWC_SH(self, r, properties):
        x = r / properties["rs"]
        beta = properties["rc"] / properties["rs"]
        rho = properties["rhos"] / ((x ** 4 + beta ** 4) ** (1/4) * (x + 1) ** 2)
        return rho

    def TNFWC_SH_hard_cutoff(self, r, properties):
        x = r / properties["rs"]
        beta = properties["rc"] / properties["rs"]
        rho_ = properties["rhos"] / ((x ** 4 + beta ** 4) ** (1/4) * (x + 1) ** 2)
        rho = np.where(r < properties["rt"], rho_, 0)
        return rho

    def TNFW(self, r, properties):
        x = r / properties["rs"]
        tau = properties["rt"] / properties["rs"]
        return properties["rhos"] / (x * (1 + x)**2) * tau**2 / (tau**2 + x**2)

    def NFW(self, r, properties):
        x = r / properties["rs"]
        return properties["rhos"] / (x * (1 + x)**2)

    def TNFW_dl(self, r, properties):
        x = r / properties["rs"]
        tau = properties["rt"] / properties["rs"]
        return 1 / (x * (1 + x)**2) * tau**2 / (tau**2 + x**2)


