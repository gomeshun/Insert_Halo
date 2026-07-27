import sys
sys.path.insert(0, "/home/stefan/Documents/uva/Thesis/Coding/main_code")

from IH_main import *
from IH_data_analysis import Insert_Halo_Data_Analysis
from IH_errors import StageException
from IH_errors import SettingsException
import numpy as np
#import matplotlib.pyplot as plt
#import time
import lenstronomy.Util.util as util
import tqdm





class IH_ABC():
    def __init__(self):
        self.priors = None
        self.accepted_parameter_values = np.array([])
        self.drawn_parameter_values = np.array([])
        self.N_simulations = 0
        self.summary_statistics = np.array([])
        self.LM = []
        


    def run_ABC(self, prior_dict, tolerance, system_class, zlens, zsource, DM_type, cone_angle_arcsec, center_x, center_y, precomputed_host, Mhost,
                obs_fluxes, flux_positions, flux_index_ref, acceptance_fraction = None, N_sim = None, geometry_shape = "DOUBLE_CONE", 
                shmfplindex = -1.9, msub_min = 7, msub_max = 11, 
                idle_run = False, host_halo_scale_factor = 0.88, redshift_scale_factor = 1.7, seed = None, randomization = True, spatial = "einasto", 
                correlated_clustering = True, LOS_normalization = 1.0, rerun_sashimi = False, nph_subhalos = 0, fix_mass_range_sh = True,
                sh_beta = 4, sh_tt_th = 1.1, sh_dz = 0.01, sh_N_herm = 20, sh_zmax = 10., sh_N_ma = 500, reset_data = False, einasto_alpha = 0.678, einasto_r2 = None, 
                sigmalogc = 0.128, gdecom_P = 8, ngaussian = 15, formation_mass_fraction = 0.5):
        
        """
        Main function to run the ABC method found in https://arxiv.org/pdf/1712.04945 using an Insert Halo dark matter model (either CDM or SIDM).
        This code runs using the following steps:


        - Determine prior distributions and DM type to test (to be determined by user in self.priors)
        - Load in (or simulate) observational/mock strong lensing data
        - (Optional) Fit a main lens model (host + galaxy) to the strong lensing data or read the main lens model data, based on the image positions and/or time delays
        - Randomly draw from the prior distributions and run Insert Halo and generate a substructure lensmodel in Lenstronomy
        - Add main lens model and substructure lens models together
        - Run Ray tracing to determine flux ratios for each image
        - Compute summary statistics and compute (Euclidean) distance between simulated and actual summary statistics
        - Accept or reject prior distribution draw, depending on the (Euclidean) distance and tolerance
        - Repeat until the distribution of accepted samples is stable to changes in the total number of samples computed and the total number of samples retained in posterior or
        N_sim is reached
        - Profit



        Current version will simply record all values drawn from the priors and their summary statistics.
        """

        # These first steps are to initialize the ABC procedure.

        # Setting user-determined priors:
        self.priors = prior_dict
        self.obs_flux_values = obs_fluxes
        self.flux_grid_positions = flux_positions
        self.obs_flux_ratios = self.compute_flux_ratios(self.obs_flux_values, flux_index_ref)
        self.acceptance_fraction = acceptance_fraction

        #===============================================================
        """The following should be done inside the main ABC loop"""
        #===============================================================

        if N_sim != None:
            for sim in tqdm.tqdm(range(N_sim), total = N_sim, desc = "Running ABC pipeline"):
                self.draw_from_priors()
                self.drawn_parameter_values = np.append(self.drawn_parameter_values, self.current_draws)

                self.IH = Insert_Halo()
                # This works assuming the names of the priors are the same as the run_insert_halo parameters:
                # (Still have to add values that should be known from data (such as zlens or cone_angle_arcsec))
                base_run = self.IH.run_insert_halo(zlens = zlens, zsource = zsource, DM_type = DM_type, cone_angle_arcsec = cone_angle_arcsec,
                                            center_x = center_x, center_y = center_y, Mhost = Mhost, 
                                            geometry_shape = geometry_shape, shmfplindex = shmfplindex, msub_min = msub_min, msub_max = msub_max, 
                                            idle_run = idle_run, host_halo_scale_factor = host_halo_scale_factor, redshift_scale_factor = redshift_scale_factor, 
                                            seed = seed, randomization = randomization, spatial = spatial, 
                                            correlated_clustering = correlated_clustering, LOS_normalization = LOS_normalization, rerun_sashimi = rerun_sashimi, 
                                            nph_subhalos = nph_subhalos, fix_mass_range_sh = fix_mass_range_sh,
                                            sh_beta = sh_beta, sh_tt_th = sh_tt_th, sh_dz = sh_dz, sh_N_herm = sh_N_herm, sh_zmax = sh_zmax, sh_N_ma = sh_N_ma, 
                                            reset_data = reset_data, sigmalogc = sigmalogc, gdecom_P = gdecom_P, ngaussian = ngaussian, 
                                            formation_mass_fraction = formation_mass_fraction, **self.current_draws)

                #print(base_run.IHB.sigma_0)
                DA = Insert_Halo_Data_Analysis(base_run)
                LM = DA.prepare_lensmodel(mass_sheet_correction = False, include_subhalos = True, include_fieldhalos = True, include_host = False, 
                                          host_lens_profile = "EPL", force_pointmasses = False, pre_model_return = True)

                #self.LM.append([LM["kwargs_substructure"]])
                system_class.setup_lensmodel(substructure_lensmodel_data = LM)
                system_class.setup_image_model()
                self.LM.append(len(system_class.image_model.LensModel.lens_model_list))
                simulated_image = system_class.image_model.image(
                                                kwargs_lens = system_class.fitted_model_kwargs["kwargs_lens"] + LM["kwargs_substructure"],
                                                kwargs_source = system_class.fitted_model_kwargs["kwargs_source"],
                                                kwargs_lens_light = system_class.fitted_model_kwargs["kwargs_lens_light"],
                                                kwargs_ps = system_class.fitted_model_kwargs["kwargs_ps"]
                                                )

                self.sim_flux_values = system_class.get_flux_image_from_sim(simulated_image, self.flux_grid_positions, num_features = 4)
                #print(self.sim_flux_values)
                self.sim_flux_ratios = self.compute_flux_ratios(self.sim_flux_values, flux_index_ref)
                self.compute_summary_statistics()
                if sim != N_sim - 1:
                    self.reset_IH()
                    self.reset_sim_fluxes()
                #print(base_run.PH.SIDM_dist.halos[0])
                system_class.reset_image_model()
                #print(self.summary_statistics[-1])
        else:
            raise NotImplementedError("Currently the run_ABC function only supports the number of simulation runs to be given by N_sim.")
        #self.acceptance_deter(tolerance = tolerance, acceptance_fraction = acceptance_fraction)


    def load_obs_data(self):
        pass

    def load_main_lens(self):
        pass

    def prepare_parameter_acceptance(self):
        if self.priors == None:
            raise StageException("Priors have not been given, cannot initialize parameter acceptance.")
        for prior in self.priors:
            self.accepted_parameter_values[prior] = []

        
    def draw_from_priors(self):
        """
        Randomly draw data points from the given priors, which should be found in self.priors. Here we assume that the self.priors variable is a dictionary.
        """
        draws = {}
        
        for prior in self.priors:
            prior_class = self.priors[prior]
            distr_type = prior_class.distr
            prior_range = prior_class.draw_range
            if distr_type == "uniform":
                draws[prior] = np.random.uniform(low = prior_range["min"], high = prior_range["max"])
            elif distr_type == "gaussian":
                draws[prior] = np.random.normal(loc = prior_range["center"], scale = prior_range["std"])
            elif distr_type == "loguniform":
                draw_uniform = np.random.uniform(low = np.log10(prior_range["min"]), high = np.log10(prior_range["max"]))
                draws[prior] = 10**draw_uniform
        self.current_draws = draws

            
    def compute_flux_ratios(self, image_fluxes, flux_index_ref):
        """
        Compute the observed or simulated flux ratios using the flux value with index flux_index_ref as the reference flux
        """

        ref_flux = image_fluxes[flux_index_ref]
        flux_ratios = []

        for flux in image_fluxes:
            flux_ratios.append(flux / ref_flux)
            
        del flux_ratios[flux_index_ref]
        return flux_ratios


    def compute_summary_statistics(self):
        summary_stat = 0

        # This shouldn't happen, but just in case this error is here to detect a missmatch between the number of observed and simulated flux ratios
        # Computing the summary statistics should therefore be done in a try except block so that data that does raise this error can be excluded.
        # This summary statistic is introduced by https://arxiv.org/pdf/1712.04945
        if len(self.obs_flux_ratios) != len(self.sim_flux_ratios):
            raise ValueError("Summary statistics cannot be computed as the number of flux ratios differs between the observed and simulated flux ratios")

        for obs_flux_ratio, sim_flux_ratio in zip(self.obs_flux_ratios, self.sim_flux_ratios):
            summary_stat += (obs_flux_ratio - sim_flux_ratio)**2
        summary_statistic = np.sqrt(summary_stat)
        self.summary_statistics = np.append(self.summary_statistics, summary_statistic)
        #return summary_statistic


    def eucl_distance(self, comp_values, new_values):
        if len(comp_values) == len(new_values) or (len(comp_values) == 1 or len(new_values) == 1):
            comp_values = np.array(comp_values)
            new_values = np.array(new_values)
            return np.sqrt((comp_values - new_values)**2)


    def acceptance_deter(self, tolerance, acception_fraction):
        if acceptance_fraction != None:    
            pass

    def reset_sim_fluxes(self):
        if hasattr(self, "sim_flux_ratios"):
            delattr(self, "sim_flux_ratios")
        if hasattr(self, "sim_flux_values"):
            delattr(self, "sim_flux_values")
            
    def reset_IH(self):
        if hasattr(self, "IH"):
            delattr(self, "IH")

    def sort_sum_stats(self, acc_fraction = None):
        ordered_indices = np.array(sorted(range(len(self.summary_statistics)), key = lambda k: self.summary_statistics[k]))
        
        if acc_fraction != None:
            N_acc = int(len(self.summary_statistics) * acc_fraction)
            accepted_indices = ordered_indices[:N_acc]
            return accepted_indices, self.summary_statistics[accepted_indices], self.drawn_parameter_values[accepted_indices]
        else:
            return ordered_indices, self.summary_statistics[ordered_indices], self.drawn_parameter_values[ordered_indices]


    def store_ABC_data(self, path, file_name, summary_statistics = None, parameter_values = None):
        if summary_statistics == None or parameter_values == None:
            if summary_statistics == None:
                summary_statistics = self.summary_statistics
            elif parameter_values == None:
                parameter_values = self.drawn_parameter_values

        assert len(summary_statistics) == len(parameter_values)
        assert type(path) == str and type(file_name) == str
    
        file_location = path + "/" + file_name
        e_test = Path(file_location + ".npz")
        if e_test.exists():
            raise NameError("A file with this name already exists, to make sure it is not overwritten, please use a different name.")
    
        
        np.savez(file_location, summary_statistics = summary_statistics, parameter_values = parameter_values, allow_pickle = True)
        
    def get_ABC_data(self, path, file_name):
        file_location = path + "/" + file_name
        return np.load(file_location, allow_pickle = True)
        
class priorc():
    def __init__(self, prior_name, distribution_type, prior_distr_values):
        self.name = prior_name
        self.distr = distribution_type

        if self.distr == "uniform":
            try:
                min_value, max_value = prior_distr_values["min"], prior_distr_values["max"]
            except KeyError:
                raise NameError(f"Prior class of distribution type {self.distr} should have a minimum and maximum value")
        elif self.distr == "gaussian":
            try:
                center_value, std = prior_distr_values["center"], prior_distr_values["std"]
            except KeyError:
                raise NameError(f"Prior class of distribution type {self.distr} should have a center and standard deviation value")
        elif self.distr == "loguniform":
            try:
                min_value, max_value = prior_distr_values["min"], prior_distr_values["max"]
            except KeyError:
                raise NameError(f"Prior class of distribution type {self.distr} should have a minimum and maximum value")            
        self.draw_range = prior_distr_values




    
