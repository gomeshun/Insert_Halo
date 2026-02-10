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





    def plot_NFW_parameters(self, PH_comp = True, save_fig = False, save_path = None, light_pdf = True):
        if PH_comp:
            cdm_ph = self.PH.run_pyhalo_CDM(n_subhalos = None, rr_background = True)

        fig = plt.figure(1)
        ax1 = plt.subplot(141)
        ax2 = plt.subplot(142)
        ax3 = plt.subplot(143)
        ax4 = plt.subplot(144)
        fig.set_size_inches(20,4)
        #plt.suptitle(rf"NFW parameter comparison for z: {zlens}", fontsize = 16)


        masses_sh = []
        rs_sh = []
        rhos_sh = []
        rt_sh = []
        
        for halo in self.PH.cdm.subhalos:
            masses_sh.append(halo.mass)
            params = halo.params_physical
        
            rs_sh.append(params["rs"])
            rhos_sh.append(params["rhos"])
            rt_sh.append(params["r_trunc_kpc"])

        ax1.scatter(rs_sh, rhos_sh, label = "Sashimi", s = 3)
        ax1.set_yscale("log")
        ax1.set_xlabel(r"$r_s$ [kpc]")
        ax1.set_ylabel(r"$\rho_s$ [$M_{\odot}$ $kpc^{-3}$]")
        ax1.set_rasterized(light_pdf)


        ax2.scatter(masses_sh, rhos_sh, label = "Sashimi", s = 3)
        ax2.set_yscale("log")
        ax2.set_xscale("log")
        ax2.set_xlabel(r"Mass [$M_{\odot}$]")
        ax2.set_ylabel(r"$\rho_s$ [$M_{\odot}$ $kpc^{-3}$]")
        ax2.set_rasterized(light_pdf)


        ax3.scatter(masses_sh, rs_sh, label = "Sashimi", s = 3)
        ax3.set_yscale("log")
        ax3.set_xscale("log")
        ax3.set_xlabel(r"Mass [$M_{\odot}$]")
        ax3.set_ylabel(r"$r_s$ [kpc]")
        ax3.set_rasterized(light_pdf)


        ax4.scatter(masses_sh, rt_sh, label = "Sashimi", s = 3)
        ax4.set_yscale("log")
        ax4.set_xscale("log")
        ax4.set_xlabel(r"Mass [$M_{\odot}$]")
        ax4.set_ylabel(r"$r_t$ [kpc]")
        ax4.set_rasterized(light_pdf)


        

        if PH_comp:
            masses_ph = []
            rs_ph = []
            rhos_ph = []
            rt_ph = []
            
            for halo in cdm_ph.subhalos:
                masses_ph.append(halo.mass)
                params_phys = halo.params_physical
                
                rs_ph.append(params_phys["rs"])
                rhos_ph.append(params_phys["rhos"])
                rt_ph.append(params_phys["r_trunc_kpc"])
                
            ax1.scatter(rs_ph, rhos_ph, label = "pyHalo", s = 3)
            ax1.legend()
            ax2.scatter(masses_ph, rhos_ph, label = "pyHalo", s = 3)
            ax2.legend() 
            ax3.scatter(masses_ph, rs_ph, label = "pyHalo", s = 3)
            ax3.legend()
            ax4.scatter(masses_ph, rt_ph, label = "pyHalo", s = 3)
            ax4.legend()

        
        plt.tight_layout()        
       
        
        if save_fig:
            if save_path == None:
                raise NameError("A path and name need to be given to save this figure, this can be done with the 'save_path' argument.")
            else:
                plt.savefig(f"{save_path}", bbox_inches = "tight")
        plt.show() 
        
        

        
        

        
        

        

        





