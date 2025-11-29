from pyHalo.PresetModels.cdm import CDM
import numpy as np
from IH_errors import StageException


class PH_run_storage():
    def __init__(self, geometry_and_massfunction):
        self.GM = geometry_and_massfunction
        self.IHB = self.GM.IHB
        self.sh = self.GM.SH.sh_settings
        self.cdm = None
        self.ran_pyhalo = False

    def run_pyhalo_CDM(self, n_subhalos = 0):
        self.IHB.global_status.set_status(f"Running pyHalo-{self.IHB.DM_type}")
        print(self.IHB.global_status.status())

        if n_subhalos == None:
            sigma_sub_corrected = self.IHB.sigma_sub
        else:
            sigma_sub_corrected = n_subhalos / (self.GM.lens_plane_area() * self.GM.integrated_subhalo_mass_function()) * self.IHB.sigma_sub

        cdm =  CDM(z_lens = self.IHB.zlens, z_source = self.IHB.zsource, 
                    sigma_sub = sigma_sub_corrected, 
                    log_mlow = self.IHB.msub_min, log_mhigh = self.IHB.msub_max, log10_sigma_sub=None,
                    log10_dNdA = None, concentration_model_subhalos='LUDLOW2016', 
                    kwargs_concentration_model_subhalos={}, concentration_model_fieldhalos='LUDLOW2016', kwargs_concentration_model_fieldhalos={},
                    truncation_model_subhalos='TRUNCATION_GALACTICUS', kwargs_truncation_model_subhalos={},
                    truncation_model_fieldhalos='TRUNCATION_RN', kwargs_truncation_model_fieldhalos={},
                    infall_redshift_model='HYBRID_INFALL', kwargs_infall_model={},
                    subhalo_spatial_distribution = self.IHB.spatial.upper(),
                    shmf_log_slope = self.IHB.shmfplindex, cone_opening_angle_arcsec = self.IHB.cone_angle_arcsec,
                    log_m_host = self.IHB.Mhost,  r_tidal=0.25, LOS_normalization = self.IHB.LOS_normalization, two_halo_contribution=self.IHB.correlated_clustering,
                    delta_power_law_index=0.0, geometry_type=self.IHB.geometry_shape.upper(), 
                    kwargs_cosmo = {"H0": self.sh.h*100, "Om0": self.sh.OmegaM, "Ob0": self.sh.OmegaB}, 
                    host_scaling_factor = self.IHB.host_halo_scale_factor,
                    redshift_scaling_factor = self.IHB.redshift_scale_factor, two_halo_Lazar_correction=self.IHB.correlated_clustering, 
                    draw_poisson=False, c_host=6.0,
                    add_globular_clusters=False, kwargs_globular_clusters=None, mass_threshold_sis=5*10**10,
                    galaxy_model='GNFW', halo_mass_profile='TNFW')
        
        self.cdm = cdm
        self.ran_pyhalo = True
        #print(sigma_sub_corrected)

        self.IHB.global_status.set_status("idle")



                    
        



















    