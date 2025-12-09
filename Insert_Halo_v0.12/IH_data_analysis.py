import numpy as np
import matplotlib.pyplot as plt




class Insert_Halo_Data_Analysis():
    def __init__(self, halo_conversion_class):

        # Easy access to all data:
        self.HC = halo_conversion_class
        self.SDG = self.HC.SDG
        self.PH = self.HC.PH
        self.cosmology = self.HC.cosmology
        self.SH = self.HC.SH
        self.IHB = self.SDG.IHB
        self.GM = self.SDG.GM


    def plot_subhalo_mass_function(self, field_halos = True, nbins = 20, return_data = False):
        min_mass = self.IHB.msub_min
        max_mass = self.IHB.msub_max

        if self.IHB.DM_type == "CDM":
            DM_realization = self.PH.cdm
        elif self.IHB.DM_type == "SIDM":
            print("The dark matter type 'SIDM' has not yet been implemented")

        subhalo_masses = [halo.mass for halo in DM_realization.subhalos]
        fieldhalo_masses = [halo.mass for halo in DM_realization.field_halos]
        
        bins = np.logspace(min_mass, max_mass, nbins)
        
        Nsub, msub = np.histogram(subhalo_masses, bins = bins)
        Nf, mf = np.histogram(fieldhalo_masses, bins = bins)

        Nsubs_mass_function = self.GM.integrated_subhalo_mass_function_massbins(np.log10(msub[0:-1]), np.log10(msub[1:])) * self.GM.lens_plane_area()

        msub_mid = (msub[1:] + msub[:-1]) / 2

        plt.loglog(msub_mid, Nsubs_mass_function, label = "Subhalo mass function", color = "black", linestyle = "dashed", zorder = 5)
        plt.loglog(msub_mid, Nsub, label = "Subhalos (Sashimi)")
        if field_halos:
            plt.loglog(msub_mid, Nf, label = "Field halos (pyHalo)")
        plt.xlabel('Halo mass '+r'$M_{\odot}$', fontsize=16)
        plt.ylabel('n(M)', fontsize=16)
        plt.legend()

        if return_data:
            if field_halos:
                return msub_mid, Nsub, Nsubs_mass_function, Nf
            else:
                return msub_mid, Nsub, Nsubs_mass_function
        

    def plot_spatial_3d(self):
        fig = plt.figure()
        fig.set_size_inches(14,12)
        ax = plt.subplot(111, projection = "3d")
        self.PH.cdm.plot(ax)











