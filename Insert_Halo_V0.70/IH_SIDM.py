
class SIDM_computations():
    def __init__(self, geometry_and_massfunction):
        self.GM = geometry_and_massfunction
        self.sh = self.GM.SH.sh_settings
        self.IHB = self.GM.IHB

        self.sigma_eff_interp = self.sh.sigma_eff_m_interpolate_analytical(sigma0_m = self.IHB.sigma_0 * self.sh.cm **2 / self.sh.gram,
                                                                           w = self.IHB.vel_scale * self.sh.km / self.sh.s)
        self.collapse_time = self.sh.t_collapse


#sigma_eff_interp = IHf.base_run.sh.sigma_eff_m_interpolate_analytical(sigma0_m = sigma_0_m, w = w)
#c_times = IHf.base_run.sh.t_collapse(sigma_eff_m = sigma_eff_array, rmax = rmax_array, Vmax = vmax_array) / sh.Gyr

