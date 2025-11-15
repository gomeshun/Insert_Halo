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



# Adding a new parameter to the run_insert_halo function requires that this gets added to at least the init function of 
# global_parameters in the Insert_halo_base class.

def run_insert_halo(zlens = 0.5, zsource = 2.0, DM_type = "CDM", sigma_0 = 20, vel_scale = 25, Mhost = 13, geometry_shape = 
                   "DOUBLE_CONE", shmfplindex = -1.9, msub_min = 6, msub_max = 10, sigma_sub = 0.025, cone_angle_arcsec = 6.0, idle_run = True, 
                   host_halo_scale_factor = 0.88, redshift_scale_factor = 1.7, seed = 1234, randomization = True, spatial = "uniform", 
                   correlated_clustering = True):

    """
    Master function for Insert_Halo. 

    ================================================================
    Currently still under active development
    Code is not yet finished and will likely contain bugs

    Current version: 0.06
        Adding geometry_and_massfunction class to compute how many
        subhalos should be generated from sashimi in the given 
        geometry
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
    """

    # Initiate the Insert Halo Base class which contains all the base information required to run Insert Halo:
    IHB = Insert_halo_base(zlens = zlens, zsource = zsource, DM_type = DM_type, sigma_0 = sigma_0,
                          vel_scale = vel_scale, Mhost = Mhost, geometry_shape = geometry_shape,
                          shmfplindex = shmfplindex, msub_min = msub_min, msub_max = msub_max, 
                           sigma_sub = sigma_sub, cone_angle_arcsec = cone_angle_arcsec, 
                           host_halo_scale_factor = host_halo_scale_factor, redshift_scale_factor = redshift_scale_factor,
                          seed = seed, randomization = randomization, spatial = spatial, correlated_clustering = correlated_clustering)

    # Initiate the Sashimi class and run it using the relevant properties defined in IHB:
    SH_run = SH_control_storage(IHB)

    if not idle_run:
        SH_run.run_sashimi()
    if idle_run:
        SH_run.init_sh_settings()

    Cosmology = cosmology(SH_run)
    GM = geometry_and_massfunction(SH_run, Cosmology)
    



    PH = PH_run_storage(GM)

    #if not idle_run:
    PH.run_pyhalo_CDM(n_subhalos = 1)

    SDG = spatial_distribution_generation(PH)


    # Run pyHalo?
    # Give pyHalo class to class that adds pyHalo and current IH class


    if not idle_run:
        dist_indices = SDG.distribution_from_weights()
        SH_run.select_data(dist_indices)

    return SDG



    
    """
    Next steps:
    
        -Define shape and size of geometry (in case of cone, make sure shape and size match with circular area at zlens that matches the cone generated
        by pyHalo)
        
        -Use subhalo mass function (like https://arxiv.org/pdf/2105.05259) to generate a statistically feasibly number for the amount of subhalos one
        would expect in the defined geometry

        -Using sashimi data and the defined number of requires subhalos start the MC procedure and pick N number of subhalos from sashimi data

        -




    Notes:
        -When running pyhalo we can let it generate only subhalos or field halos by setting either LOS_normalization = 0.0 (no field halos) or sigma_sub = 0.0 (no subhalos) This will be useful for performance when everything works later.
    """


    


