import os
import sys
file_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, file_path)
from Insert_halo_base import *
from sashimi_si import *
from IH_errors import *

class SH_control_storage(Insert_halo_base):
    """
    Class which stores all of the sashimi data and outputs as well as several functions to control and call sashimi.
    """
    def __init__(self, Insert_halo_base):
        self.IHB = Insert_halo_base


        #Sashimi storage values:
        self.ma200 = None
        self.z_acc = None
        self.rsCDM_acc = None
        self.rhosCDM_acc = None
        self.rmaxCDM_acc = None
        self.VmaxCDM_acc = None
        self.rsSIDM_acc = None
        self.rhosSIDM_acc = None
        self.rcSIDM_acc = None
        self.rmaxSIDM_acc = None
        self.VmaxSIDM_acc = None
        self.m_z0 = None
        self.rsCDM_z0 = None
        self.rhosCDM_z0 = None
        self.rmaxCDM_z0 = None
        self.VmaxCDM_z0 = None
        self.rsSIDM_z0 = None
        self.rhosSIDM_z0 = None
        self.rcSIDM_z0 = None
        self.rmaxSIDM_z0 = None
        self.VmaxSIDM_z0 = None
        self.ctCDM_z0 = None
        self.tt_ratio = None
        self.weightCDM = None
        self.weightSIDM = None
        self.surviveCDM = None
        self.surviveSIDM = None

        self.sh_settings = None

        self.ran_sashimi = False



        # Selected sashimi storage values (after generating a random distribution of subhalos from the geometry and massfunction)
        self.S_ma200 = None
        self.S_z_acc = None
        self.S_rsCDM_acc = None
        self.S_rhosCDM_acc = None
        self.S_rmaxCDM_acc = None
        self.S_VmaxCDM_acc = None
        self.S_rsSIDM_acc = None
        self.S_rhosSIDM_acc = None
        self.S_rcSIDM_acc = None
        self.S_rmaxSIDM_acc = None
        self.S_VmaxSIDM_acc = None
        self.S_m_z0 = None
        self.S_rsCDM_z0 = None
        self.S_rhosCDM_z0 = None
        self.S_rmaxCDM_z0 = None
        self.S_VmaxCDM_z0 = None
        self.S_rsSIDM_z0 = None
        self.S_rhosSIDM_z0 = None
        self.S_rcSIDM_z0 = None
        self.S_rmaxSIDM_z0 = None
        self.S_VmaxSIDM_z0 = None
        self.S_ctCDM_z0 = None
        self.S_tt_ratio = None
        self.S_weightCDM = None
        self.S_weightSIDM = None
        self.S_surviveCDM = None
        self.S_surviveSIDM = None

    def init_sh_settings(self):
        """
        Initializes Sashimi and returns a Sashimi-SIDM object without running it. This is used for development purposes.
        """
        if self.IHB.global_status.status() != "idle":
            raise StatusException(r"Global program status should be idle, but is set to {self.IHB.global_status.status()} in init_sh_settings().")

        #Setting global status:
        self.IHB.global_status.set_status(f"Initializing Sashimi-{self.IHB.DM_type}")
        print(self.IHB.global_status.status())

        sh = subhalo_properties(sigma0_m = self.IHB.sigma_0, w = self.IHB.vel_scale, beta = 4, tt_th = 1.1)
        self.sh_settings = sh
        self.IHB.global_status.set_status("idle")
        
        
    def run_sashimi(self):
        """
        Runs Sashimi-SIDM using the properties defined in the global_parameters class which is a base class for Insert_halo_base.


        List of things to be fixed/checked:

        -sigma0_m is now given from IHB.sigma_0, but sigma_0 is not defined with DM particle mass
        -In subhalo_properties_calc the M0 (or host mass) mass is defined to be M_200, but the one given is IHB.Mhost which is assumed
        to be the total host halo mass (just like it is in pyHalo). Also in pyHalo this should be given in a base of log10
        This might be the same for pyHalo where the host halo mass is defined as M200

        
        -Control of the other (default) parameters of sashimi should be added such that sashimi can be called fully with the 
        Insert_halo_base class
        
        """
        #Be careful with passing sigma0_m as self.IHB.sigma_0 here as sigma_0 is not defined with the DM particle mass unlike sigma0_m

        #Checking global status
        if self.IHB.global_status.status() != "idle":
            raise StatusException(r"Global program status should be idle, but is set to {self.IHB.global_status.status()} in run_sashimi().")

        #Setting global status:
        self.IHB.global_status.set_status(f"Running Sashimi-{self.IHB.DM_type}")
        print(self.IHB.global_status.status())


        sh = subhalo_properties(sigma0_m = self.IHB.sigma_0, w = self.IHB.vel_scale, beta = 4, tt_th = 1.1)
        
        ma200, z_acc, rsCDM_acc, rhosCDM_acc, rmaxCDM_acc, VmaxCDM_acc, rsSIDM_acc, rhosSIDM_acc, rcSIDM_acc, rmaxSIDM_acc, VmaxSIDM_acc, m_z0, rsCDM_z0,   rhosCDM_z0, rmaxCDM_z0, VmaxCDM_z0, rsSIDM_z0, rhosSIDM_z0, rcSIDM_z0, rmaxSIDM_z0, VmaxSIDM_z0, ctCDM_z0, tt_ratio, weightCDM, weightSIDM, surviveCDM, surviveSIDM = sh.subhalo_properties_calc(M0 = 10 ** self.IHB.Mhost, redshift = self.IHB.zlens, M0_at_redshift = True, dz = 0.01, N_herm = 20, zmax = 5.,                            logmamin = self.IHB.msub_min, logmamax = self.IHB.msub_max, N_ma = 500)

        self.ma200 = ma200
        self.z_acc = z_acc
        self.rsCDM_acc = rsCDM_acc
        self.rhosCDM_acc = rhosCDM_acc
        self.rmaxCDM_acc = rmaxCDM_acc
        self.VmaxCDM_acc = VmaxCDM_acc
        self.rsSIDM_acc = rsSIDM_acc
        self.rhosSIDM_acc = rhosSIDM_acc
        self.rcSIDM_acc = rcSIDM_acc
        self.rmaxSIDM_acc = rmaxSIDM_acc
        self.VmaxSIDM_acc = VmaxSIDM_acc
        self.m_z0 = m_z0
        self.rsCDM_z0 = rsCDM_z0
        self.rhosCDM_z0 = rhosCDM_z0
        self.rmaxCDM_z0 = rmaxCDM_z0
        self.VmaxCDM_z0 = VmaxCDM_z0
        self.rsSIDM_z0 = rsSIDM_z0
        self.rhosSIDM_z0 = rhosSIDM_z0
        self.rcSIDM_z0 = rcSIDM_z0
        self.rmaxSIDM_z0 = rmaxSIDM_z0
        self.VmaxSIDM_z0 = VmaxSIDM_z0
        self.ctCDM_z0 = ctCDM_z0
        self.tt_ratio = tt_ratio
        self.weightCDM = weightCDM
        self.weightSIDM = weightSIDM
        self.surviveCDM = surviveCDM
        self.surviveSIDM = surviveSIDM

        self.sh_settings = sh

        self.ran_sashimi = True
        
        self.IHB.global_status.set_status("idle")
 

    def renormalize_weights(self, N_new_total):

        if not self.ran_sashimi:
            print("Warning, Sashimi has not been run yet, returning None")
            return None
        elif self.ran_sashimi:
            if self.IHB.DM_type == "CDM":
                w_new = self.weightCDM * N_new_total / np.sum(self.weightCDM)
            elif self.IHB.DM_type == "SIDM":
                w_new = self.weightSIDM * N_new_total / np.sum(self.weightSIDM)
            return w_new


    def select_data(self, indices):
        self.S_ma200 = self.ma200[indices]
        self.S_z_acc = self.z_acc[indices]
        self.S_rsCDM_acc = self.rsCDM_acc[indices]
        self.S_rhosCDM_acc = self.rhosCDM_acc[indices]
        self.S_rmaxCDM_acc = self.rmaxCDM_acc[indices]
        self.S_VmaxCDM_acc = self.VmaxCDM_acc[indices]
        self.S_rsSIDM_acc = self.rsSIDM_acc[indices]
        self.S_rhosSIDM_acc = self.rhosSIDM_acc[indices]
        self.S_rcSIDM_acc = self.rcSIDM_acc[indices]
        self.S_rmaxSIDM_acc = self.rmaxSIDM_acc[indices]
        self.S_VmaxSIDM_acc = self.VmaxSIDM_acc[indices]
        self.S_m_z0 = self.m_z0[indices]
        self.S_rsCDM_z0 = self.rsCDM_z0[indices]
        self.S_rhosCDM_z0 = self.rhosCDM_z0[indices]
        self.S_rmaxCDM_z0 = self.rmaxCDM_z0[indices]
        self.S_VmaxCDM_z0 = self.VmaxCDM_z0[indices]
        self.S_rsSIDM_z0 = self.rsSIDM_z0[indices]
        self.S_rhosSIDM_z0 = self.rhosSIDM_z0[indices]
        self.S_rcSIDM_z0 = self.rcSIDM_z0[indices]
        self.S_rmaxSIDM_z0 = self.rmaxSIDM_z0[indices]
        self.S_VmaxSIDM_z0 = self.VmaxSIDM_z0[indices]
        self.S_ctCDM_z0 = self.ctCDM_z0[indices]
        self.S_tt_ratio = self.tt_ratio[indices]
        self.S_weightCDM = self.weightCDM[indices]
        self.S_weightSIDM = self.weightSIDM[indices]
        self.S_surviveCDM = self.surviveCDM[indices]
        self.S_surviveSIDM = self.surviveSIDM[indices]

"""
For CDM:
        ma200:    Mass m_{200} at accretion.
        z_acc:    Redshift at accretion.
        rs_acc:   Scale radius r_s at accretion.
        rhos_acc: Characteristic density \rho_s at accretion.
        m_z0:     Mass up to tidal truncation radius at a given redshift.
        rs_z0:    Scale radius r_s at a given redshift.
        rhos_z0:  Characteristic density \rho_s at a given redshift.
        ct_z0:    Tidal truncation radius in units of r_s at a given redshift.
        weight:   Effective number of subhalos that are characterized by the same set of the parameters above.
        survive:  If that subhalo survive against tidal disruption or not.
"""








