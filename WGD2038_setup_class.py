import numpy as np
import copy
from scipy import ndimage



class WGD2038():
    def __init__(self):
        self.data = None
        self.model = None
        self.kwargs_model = None
        self.fitted_model_kwargs = None
        self.usable_lensmodel = None
        self.usable_sourcemodel = None
        self.usable_lenslight_model = None


        # Imports from Samana and Lenstronomy:
        
        from samana.Data.wgd2038 import WGD2038_HST
        from samana.Model.wgd2038_model import WGD2038ModelEPLM3M4Shear
        from lenstronomy.Workflow.fitting_sequence import FittingSequence
        from lenstronomy.Data.coord_transforms import Coordinates
        from lenstronomy.Plots.model_plot import ModelPlot
        from lenstronomy.LensModel.lens_model import LensModel
        from lenstronomy.LightModel.light_model import LightModel
        from lenstronomy.Data.imaging_data import ImageData
        from lenstronomy.Data.psf import PSF
        from lenstronomy.ImSim.image_model import ImageModel
        from lenstronomy.PointSource.point_source import PointSource




        # Insert_Halo imports:
        from IH_errors import StageException, SettingsException
        
        self.WGD2038_HST = WGD2038_HST
        self.WGD2038ModelEPLM3M4Shear = WGD2038ModelEPLM3M4Shear
        self.FittingSequence = FittingSequence
        self.Coordinates = Coordinates
        self.ModelPlot = ModelPlot
        self.LensModel = LensModel
        self.LightModel = LightModel
        self.ImageData = ImageData
        self.PSF = PSF
        self.ImageModel = ImageModel
        self.PointSource = PointSource
        self.StageException = StageException
        self.SettingsException = SettingsException

        
    def general_setup(self):
        self.data = self.WGD2038_HST()
        self.model = self.WGD2038ModelEPLM3M4Shear(self.data, shapelets_order=16)



        
    def perform_PSO(self, cosmology):
        self.data = self.WGD2038_HST()
        self.model = self.WGD2038ModelEPLM3M4Shear(self.data, shapelets_order=16)
        kwargs_lens_macro_init = None
        kwargs_model, _, _, _, _ = self.model.setup_kwargs_model(kwargs_lens_macro_init=kwargs_lens_macro_init, astropy_cosmo = cosmology)

        self.kwargs_model = kwargs_model
        kwargs_params = self.model.kwargs_params(kwargs_lens_macro_init=kwargs_lens_macro_init)
        fitting_sequence = self.FittingSequence(self.data.kwargs_data_joint, self.kwargs_model, self.model.kwargs_constraints, 
                                              self.model.kwargs_likelihood, kwargs_params)
        nthreads = 10

        fitting_kwargs_list = [
                ['update_settings', {'lens_add_fixed': [[0, ['a1_a', 'a4_a', 'a3_a'], [0.0, 0.0, 0.0]]]}],
                ['PSO', {'sigma_scale': 1., 'n_particles': 10, 'n_iterations': 100, 'threadCount': nthreads}],]

        chain_list = fitting_sequence.fit_sequence(fitting_kwargs_list)
        self.fitted_model_kwargs = fitting_sequence.best_fit()

        self.image_width_arcsec = self.kwargs_model["kwargs_lens_eqn_solver"]["search_window"]
        self.x_center = self.kwargs_model["kwargs_lens_eqn_solver"]["x_center"]
        self.y_center = self.kwargs_model["kwargs_lens_eqn_solver"]["y_center"]


    def setup_data_modelplot(self):
        if not hasattr(self, "multi_band_list"):
            self.multi_band_list = self.data.kwargs_data_joint["multi_band_list"]
        self.modelplot_ns = self.ModelPlot(self.multi_band_list, self.kwargs_model, self.fitted_model_kwargs, arrow_size = 0.02,
                                          cmap_string = "bone", image_likelihood_mask_list = [self.data.likelihood_mask], fast_caustic = True)

    def setup_sub_modelplot(self, substructure_lensmodel_data):
        if not hasattr(self, "multi_band_list"):
            self.multi_band_list = self.data.kwargs_data_joint["multi_band_list"]

        fitted_model_kwargs_c = copy.deepcopy(self.fitted_model_kwargs)
        kwargs_model_c = copy.deepcopy(self.kwargs_model)

        fitted_model_kwargs_c["kwargs_lens"] = fitted_model_kwargs_c["kwargs_lens"] + substructure_lensmodel_data["kwargs_substructure"]
        kwargs_model_c["lens_model_list"] = kwargs_model_c["lens_model_list"] + substructure_lensmodel_data["lens_model_list"]
        kwargs_model_c["lens_redshift_list"] = kwargs_model_c["lens_redshift_list"] + substructure_lensmodel_data["lens_redshift_array"]

        self.modelplot_s = self.ModelPlot(self.multi_band_list, kwargs_model_c, fitted_model_kwargs_c, arrow_size = 0.02,
                                          cmap_string = "bone", image_likelihood_mask_list = [self.data.likelihood_mask], fast_caustic = True)    

    def setup_lensmodel(self, substructure_lensmodel_data = None, data_only = False):
        if data_only:
            lens_model_class = self.LensModel(lens_model_list = self.kwargs_model["lens_model_list"],
                                              lens_redshift_list = self.kwargs_model["lens_redshift_list"],
                                              multi_plane = self.kwargs_model["multi_plane"],
                                              z_source = self.kwargs_model["z_source"])
        elif not data_only:
            if substructure_lensmodel_data == None:
                raise ValueError("To add the substructure data to the lensmodel, substructure_lensmodel_data needs to be given")
            lens_model_class = self.LensModel(lens_model_list = self.kwargs_model["lens_model_list"] + substructure_lensmodel_data["lens_model_list"],
                                              lens_redshift_list = self.kwargs_model["lens_redshift_list"] + substructure_lensmodel_data["lens_redshift_array"],
                                              multi_plane = self.kwargs_model["multi_plane"],
                                              z_source = self.kwargs_model["z_source"])     
        self.usable_lensmodel = lens_model_class
        
            
    def setup_lightmodels(self, add_lens_light = False, add_shapelets = False):
        if add_shapelets:
            self.usable_sourcemodel = self.LightModel(light_model_list = self.kwargs_model["source_light_model_list"])
        elif not add_shapelets:
            print("Removing shapelets from the source model is currently not implemented.")
        if add_lens_light:
            self.usable_lenslight_model = self.LightModel(light_model_list = self.kwargs_model["lens_light_model_list"])
            

    def setup_image_model(self):
        if self.usable_sourcemodel == None:
            print("Warning, no source model was given")
        kwargs_data, kwargs_psf, kwargs_numerics = self.data.kwargs_data_joint["multi_band_list"][0]
        data_class = self.ImageData(**kwargs_data)
        psf_class = self.PSF(**kwargs_psf)

        point_source_class = self.PointSource(point_source_type_list = self.kwargs_model["point_source_model_list"], fixed_magnification_list=[False])

        self.image_model = self.ImageModel(
                                        data_class = data_class,
                                        psf_class = psf_class,
                                        lens_model_class = self.usable_lensmodel,
                                        source_model_class = self.usable_sourcemodel,
                                        lens_light_model_class = self.usable_lenslight_model,
                                        kwargs_numerics = kwargs_numerics,
                                        point_source_class = point_source_class
        )
    def reset_image_model(self):
        if hasattr(self, "image_model"):
            delattr(self, "image_model")

    def get_flux_image_from_data(self, flux_threshold_log = 0.2, data_clean_limit = 1e-10, label_galaxy = 3):
        if not hasattr(self, "multi_band_list"):
            raise self.StageException("Cannot find multi_band_list data for computing flux values, please run setup_data_modelplot() first.")
        raw_data = self.multi_band_list[0][0]["image_data"]

        raw_data_c = copy.deepcopy(raw_data)
        raw_data_c = np.where(raw_data_c <= 0, data_clean_limit, raw_data_c)
        
        mask = raw_data_c >= 10**flux_threshold_log
        labeled_array, num_features = ndimage.label(mask)

        if num_features != 4 and num_features != 5:
            raise self.SettingsException("Problem with finding the data images, likely caused by selected flux_threshold.")


        fluxes = []

        if num_features == 4:
            for num in range(num_features):
                fluxes.append(np.sum(raw_data_c[labeled_array == num + 1]))
                
            for label in range(num_features):
                    raw_data_c[labeled_array == label + 1] = np.nan   

                
        elif num_features == 5:
            labeled_array = np.where(labeled_array == label_galaxy, 0, labeled_array)
            labeled_array = np.where(labeled_array == 5, label_galaxy, labeled_array)            

            for num in range(num_features - 1):
                #if num + 1 != label_galaxy:
                fluxes.append(np.sum(raw_data_c[labeled_array == num + 1]))   

            for label in range(num_features - 1):
                #if label + 1 != label_galaxy:
                raw_data_c[labeled_array == label + 1] = np.nan      




        return fluxes, labeled_array, raw_data_c, raw_data




    def get_flux_image_from_sim(self, simulated_image, labeled_array, num_features = 4):
        assert np.max(labeled_array) == num_features

        sim_image_c = copy.deepcopy(simulated_image)

        fluxes_s = []
        
        for num in range(num_features):
            fluxes_s.append(np.sum(simulated_image[labeled_array == num + 1]))

        return fluxes_s





