import os
import sys

file_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, file_path)

from IH_run_sashimi import *
from Insert_halo_base import *
from Geometry_and_mass_function import geometry_and_massfunction
from IH_cosmology import IH_cosmology as cosmology
from IH_spatial_calc import *
from IH_run_pyhalo import *
from IH_halo_conversion import *
from IH_errors import SettingsException
from IH_errors import StageException
from IH_SIDM import SIDM_computations



# Adding a new parameter to the run_insert_halo function requires that this gets added to at least the init function of 
# global_parameters in the Insert_halo_base class as well as the rerun_insert_halo function!

class Insert_Halo():
    """
    Main Insert Halo class which contains the main function to run the program. For more information see help(Insert_Halo.run_insert_halo()).
    """
    def __init__(self):
        self.ran_Insert_Halo = False
        self.base_run = None
        self.subsequent_run = None
        

    def run_insert_halo(self, zlens = 0.5, zsource = 2.0, DM_type = "CDM", sigma_0 = 20, vel_scale = 25, Mhost = 13, geometry_shape = 
                       "DOUBLE_CONE", shmfplindex = -1.9, msub_min = 6, msub_max = 10, sigma_sub = 0.025, cone_angle_arcsec = 6.0, idle_run = True, 
                       host_halo_scale_factor = 0.88, redshift_scale_factor = 1.7, seed = None, randomization = True, spatial = "uniform", 
                       correlated_clustering = True, LOS_normalization = 1.0, rerun_sashimi = False, nph_subhalos = 0, fix_mass_range_sh = False,
                       sh_beta = 4, sh_tt_th = 1.1, sh_dz = 0.01, sh_N_herm = 20, sh_zmax = 10., sh_N_ma = 500, center_x = 0, center_y = 0,
                       reset_data = False, einasto_alpha = 0.678, einasto_r2 = 199, sigmalogc = 0.128):
    
        """
        Master function for Insert_Halo. 
    
        ================================================================
        Currently still under active development
        Code is not yet finished and will likely contain bugs
    
        Current version: 0.30
            -Added code to compute the NFW parameters of the host halo
            in a host halo class
            -Added code te add the host halo lensmodel in Lenstronomy
            based on the NFW parameters of the host for both an NFW
            lensmodel and EPL lensmodel
            -Added code to automatically generate a Lenstronomy lensmodel
            containing the desired halos
        ================================================================
        
        This function goes through the entire insert_halo procedure and generates subhalos using Sashimi. These are then given
        spatial positions alongside adjusted density profiles and cross-section models to be given to pyHalo. Then runs pyHalo
        to generate field halos (or line-of-sight halos) so that a full distribution of halos exists that can be used for
        gravitational lensing analysis using Lenstronomy.
    
    
        ===============
        Parameters
        ===============
        zlens: Redshift of the main lens plane, this is the redshift at which all subhalos are located
        zsource: Redshift of the source light behind the lensing plane(s)
        DM_type: Type of dark matter we consider (currenly only CDM is allowed, but SIDM will follow)
        sigma_0: Cross-section amplitude in [cm^2 g^-1]
        vel_scale: Velocity scale for the cross-section model [km/s]
        Mhost: log10 of the mass of host halo (M_200) in [M_sun]
        return_results: Determines what run_insert_halo returns (options for now: 'sashimi' and 'geometry')
        geometry_shape: Sets the shape of the geometry under consideration (options: 'DOUBLE_CONE', 'CYLINDER', 'CONE')
        shmfplindex: power-law index of the subhalo mass function (alpha in https://arxiv.org/pdf/1908.06983)
        msub_min: log10 of the minimal subhalo mass [log10(M_sun)]
        msub_max: log10 of the maximum subhalo mass [log10(M_sun)]
        sigma_sub: amplitude of the subhalo mass function at a mass of 10^8 solar masses [# halos / kpc^2]
        cone_angle_arcsec: Opening angle of the geometry [arcsec]
        idle_run: If True only generates an insert_halo object without actually running any code 
        host_halo_scale_factor: scale factor for the host halo (k1 in https://arxiv.org/pdf/1908.06983)
        redshift_scale_factor: scale factor for redshift (k2 in https://arxiv.org/pdf/1908.06983)
        seed: random seed used by the randomization parts of the code for result recreation. If set to None, uses standard numpy randomization
        randomization: boolean parameter which determines whether any randomization takes place (where relevant)
        spatial: str indicating the type of subhalo distribution (current option: 'uniform')
        correlated_clustering: boolean indicating whether or not to take into account correlated clustering (2-halo term) in the LOS halo mass function
        LOS_normalization: Line Of Sight (LOS) halo mass function normalization (set this to 0 to exclude all LOS halos)
        rerun_sashimi: Experimental options, currently not working properly
        nph_subhalos: [ph] Number of subhalos to be generated by pyHalo, default is 0.
        fix_mass_range_sh: [sh] boolean on whether or not the data generated by sashimi should be cut off outside the given subhalo masses msub_min and msub_max
        sh_beta: [sh] Sashimi input representing the power index of the density profile, default is 4
        sh_tt_th: [sh] Sashimi input indicating the maximum t/t_c that is simulated
        sh_dz: [sh] Sashimi input indicating redshift grid for halo accretion
        sh_N_herm: [sh] Sashimi input number of grid for the Gauss-Hermite quadrature for the integration over concentration
        sh_zmax: [sh] Sashimi input indicating the maximum input for the redshift where accretion calculation starts
        sh_N_ma: [sh] Sashimi input indicating the number of gridpoint for the subhalo mass at accretion in M_200
        center_x: x-coordinate offset of the host halo and thus the subhalo spatial distribution center in units of [arcsec]
        center_y: y-coordinate offset of the host halo and thus the subhalo spatial distribution center in units of [arcsec]
        reset_data: Boolean value indicating whether to replace the base_run data when calling run_insert_halo consecutively.
    
    
        ===============
        To add/adjust
        ===============
        -Add a way to have the cone_angle_arcsec be 3*R_Ein, where R_Ein is the Einstein radius for the given lens such as in https://arxiv.org/pdf/1901.11031
        -Add a plotting class to easily generate certain plots that are used often (halo mass function, halo spatial distribution, convergence plots)
        -Add code structure to allow the program to run multiple times without having to run all parts so it reuses some data (like a Sashimi run)
        -

        """

        # Check for resetting computed data:
        if self.ran_Insert_Halo and not reset_data:
            raise SettingsException("Insert Halo was already run, running it again will replace the stored data. If the stored data is supposed to be "
                                   "replaced, set the reset_data to True, otherwise use the rerun_insert_halo function.")
        
        
    
        # Initiate the Insert Halo Base class which contains all the base information required to run Insert Halo:
        IHB = Insert_halo_base(zlens = zlens, zsource = zsource, DM_type = DM_type, sigma_0 = sigma_0,
                              vel_scale = vel_scale, Mhost = Mhost, geometry_shape = geometry_shape,
                              shmfplindex = shmfplindex, msub_min = msub_min, msub_max = msub_max, 
                              sigma_sub = sigma_sub, cone_angle_arcsec = cone_angle_arcsec, 
                              host_halo_scale_factor = host_halo_scale_factor, redshift_scale_factor = redshift_scale_factor,
                              seed = seed, randomization = randomization, spatial = spatial, correlated_clustering = correlated_clustering,
                              LOS_normalization = LOS_normalization, fix_mass_range_sh = fix_mass_range_sh, 
                              sh_beta = sh_beta, sh_tt_th = sh_tt_th, sh_dz = sh_dz, sh_N_herm = sh_N_herm, sh_zmax = sh_zmax, sh_N_ma = sh_N_ma,
                              center_x = center_x, center_y = center_y, einasto_alpha = einasto_alpha, einasto_r2 = einasto_r2, sigmalogc = sigmalogc)
    
        # Initiate the Sashimi class and run it using the relevant properties defined in IHB:
        SH_run = SH_control_storage(IHB)
    
        if not idle_run and not SH_run.ran_sashimi:
            SH_run.run_sashimi()
        elif rerun_sashimi:
            SH_run.run_sashimi()
        if idle_run:
            SH_run.init_sh_settings()
    
        #if SH_run.ran_sashimi:
        #    SH_run.IHB.msub_min = np.log10(np.min(SH_run.m_z0))
        #    SH_run.IHB.msub_max = np.log10(np.max(SH_run.m_z0))
    
        Cosmology = cosmology(SH_run)
        GM = geometry_and_massfunction(SH_run, Cosmology)

        if IHB.DM_type == "SIDM":
            SIDM_calc = SIDM_computations(GM)
            PH = PH_run_storage(SIDM_calc)
            PH.run_pyhalo(n_subhalos = 0)
            PH.run_pyhalo_alternative()
    
        elif DM_type == "CDM":
            PH = PH_run_storage(GM)
            PH.run_pyhalo(n_subhalos = nph_subhalos)            

        SDG = spatial_distribution_generation(PH)
        HC = HaloConversion(SDG)
    
        if not idle_run:
            dist_indices = HC.SDG.distribution_from_weights()
            HC.SH.select_data(dist_indices)
            #HC.SDG.generate_uniform_distribution(len(dist_indices), x_offset = center_x, y_offset = center_y)
            x, y = HC.SDG.generate_spatial_distribution()
            HC.generate_new_halo_classes() 

            if IHB.DM_type == "CDM":
                HC.insert_CDM_halos_corr()
            elif IHB.DM_type == "SIDM":
                HC.insert_SIDM_halos()

        self.base_run = HC
        self.ran_Insert_Halo = True
        return self.base_run


    

    def rerun_insert_halo(self, zlens = 0.5, zsource = 2.0, DM_type = "CDM", sigma_0 = 20, vel_scale = 25, Mhost = 13, geometry_shape = 
                       "DOUBLE_CONE", shmfplindex = -1.9, msub_min = 6, msub_max = 10, sigma_sub = 0.025, cone_angle_arcsec = 6.0, idle_run = True, 
                       host_halo_scale_factor = 0.88, redshift_scale_factor = 1.7, seed = None, randomization = True, spatial = "uniform", 
                       correlated_clustering = True, LOS_normalization = 1.0, nph_subhalos = 0, fix_mass_range_sh = False,
                       sh_beta = 4, sh_tt_th = 1.1, sh_dz = 0.01, sh_N_herm = 20, sh_zmax = 10., sh_N_ma = 500, center_x = 0, center_y = 0,
                       replace_data = True, keep_subhalo_data = False, einasto_alpha = 0.678, einasto_r2 = 199, sigmalogc = 0.128):

        """
        Effectively the same function as run_insert_halo, but this function potentially reuses some of the data already generated by 
        running the run_insert_halo function, causing subsequent runs to take less time to compute. Typically this means that when the 
        right settings are the same, Sashimi is not run again, but the Sashimi data of a previous run is used.

        The parameters of this function are the same as the run_insert_halo function with the exception of a few settings meant specifically for 
        this rerun function, such as:

        replace_data: boolean value indicating whether the rerun data needs to replace the original data in the main data storage self.base_run
        keep_subhalo_data: boolean value indicating whether the subhalo data from the original run should be used for this one rather than
            generating a new sub-distribution from the Sashimi data like usual (This option is currently not functional)
        """


        if not self.ran_Insert_Halo:
            raise StageException("Insert Halo has not yet been run, so it cannot be rerun.")

        # Comparing the settings of the first run with that of this run:
        rerun_sashimi = False
        p_run = self.base_run
        IHB = p_run.SDG.IHB
        truth_array = np.array([IHB.zlens != zlens, IHB.DM_type != DM_type, IHB.sigma_0 != sigma_0, IHB.vel_scale != vel_scale, IHB.Mhost != Mhost,
                                IHB.msub_min != msub_min, IHB.msub_max != msub_max, IHB.fix_mass_range_sh != fix_mass_range_sh,
                                IHB.sh_beta != sh_beta, IHB.sh_tt_th != sh_tt_th, IHB.sh_dz != sh_dz, IHB.sh_N_herm != sh_N_herm,
                                IHB.sh_N_ma != sh_N_ma, IHB.sigmalogc != sigmalogc])
        if np.any(truth_array):
            rerun_sashimi = True    
        
        IHB = Insert_halo_base(zlens = zlens, zsource = zsource, DM_type = DM_type, sigma_0 = sigma_0,
                              vel_scale = vel_scale, Mhost = Mhost, geometry_shape = geometry_shape,
                              shmfplindex = shmfplindex, msub_min = msub_min, msub_max = msub_max, 
                              sigma_sub = sigma_sub, cone_angle_arcsec = cone_angle_arcsec, 
                              host_halo_scale_factor = host_halo_scale_factor, redshift_scale_factor = redshift_scale_factor,
                              seed = seed, randomization = randomization, spatial = spatial, correlated_clustering = correlated_clustering,
                              LOS_normalization = LOS_normalization, fix_mass_range_sh = fix_mass_range_sh, 
                              sh_beta = sh_beta, sh_tt_th = sh_tt_th, sh_dz = sh_dz, sh_N_herm = sh_N_herm, sh_zmax = sh_zmax, sh_N_ma = sh_N_ma,
                              center_x = center_x, center_y = center_y, einasto_alpha = einasto_alpha, einasto_r2 = einasto_r2, sigmalogc = sigmalogc)
    
        # Initiate the Sashimi class and run it using the relevant properties defined in IHB:
        SH_run = SH_control_storage(IHB)
    
        if not idle_run and rerun_sashimi:
            SH_run.run_sashimi()
        elif idle_run and not rerun_sashimi:
            SH_run.init_sh_settings()
        elif not idle_run and not rerun_sashimi:
            # To Do: Change this so that the changes to the IHB parameters can be different for the rerun
            SH_run = p_run.SH
            SH_run.IHB = IHB


        
        #if SH_run.ran_sashimi:
        #    SH_run.IHB.msub_min = np.log10(np.min(SH_run.m_z0))
        #    SH_run.IHB.msub_max = np.log10(np.max(SH_run.m_z0))
    
        Cosmology = cosmology(SH_run)
        GM = geometry_and_massfunction(SH_run, Cosmology)
        
        PH = PH_run_storage(GM)
    
    
        PH.run_pyhalo(n_subhalos = nph_subhalos)
    
        SDG = spatial_distribution_generation(PH)
        HC = HaloConversion(SDG)
    
        if not idle_run:
            dist_indices = HC.SDG.distribution_from_weights()
            HC.SH.select_data(dist_indices)
            #HC.SDG.generate_uniform_distribution(len(dist_indices), x_offset = center_x, y_offset = center_y)
            x, y = HC.SDG.generate_spatial_distribution()
            HC.generate_new_halo_classes() 

            if IHB.DM_type == "CDM":
                HC.insert_CDM_halos_corr()
            elif IHB.DM_type == "SIDM":
                HC.insert_SIDM_halos()

        if replace_data:
            self.base_run = HC
            self.ran_Insert_Halo = True
            return self.base_run
        elif not replace_data:
            self.subsequent_run = HC
            return self.subsequent_run


        



