import numpy as np

class Host_Halo():
    def __init__(self, z, mass, geometry_and_mass_function):
        self._z = z
        self._mass = mass # This should be log10(M_200)
        self.profile = "NFW"

        self.GM = geometry_and_mass_function
        self.sh = self.GM.SH.sh_settings

        

    @property
    def properties(self):
        """
        Returns the NFW_properties for the host halo using the method described in https://arxiv.org/pdf/1903.11427 and using the concentration-
        mass fitting relation from https://arxiv.org/pdf/1502.00391
        This function assumes the average concentration for the host halo (so we disregard scatter in the concentration log-normal distribution)


        returns: c200, rs, rhos, r200
        """

        if hasattr(self, "_properties"):
            simulated_radius = self.GM.kpc_per_arcsec() * self.GM.IHB.cone_angle_arcsec / 2
            if simulated_radius > self._properties["r200"]:
                print(f"Warning! The physical size of the host halo indicated by r_200 is {self._properties["r200"]:.2f} kpc which is smaller than the simulated region which is: {simulated_radius:.2f} kpc.")
            return self._properties
        else:
            c200 = self.c200()
            rs = self.get_r200_from_m200() / c200
            rhos = self.mass / (4 * np.pi * rs ** 3 * self.f_NFW(c200))
            properties = {"c200": c200, "rs": rs, "rhos": rhos, "r200": rs * c200}
            self._properties = properties

            simulated_radius = self.GM.kpc_per_arcsec() * self.GM.IHB.cone_angle_arcsec / 2
            if simulated_radius > self._properties["r200"]:
                print(f"warning! The physical size of the host halo indicated by r_200 is {self._properties["r200"]:.2f} kpc which is smaller than the simulated region which is: {simulated_radius:.2f} kpc.")
            return self._properties
        


    @property
    def logmass(self):
        return self._mass

    @property
    def mass(self):
        return 10 ** self._mass

    @property
    def z(self):
        return self._z

    def get_r200_from_m200(self):
        Delta_200 = 200
        rho_crit_z = self.sh.rhocrit(self.z) / self.sh.Msun * self.sh.kpc**3
        r_200 = (3 * self.mass / (4 * np.pi * Delta_200 * rho_crit_z)) ** (1/3)
        return r_200

    def c200(self):
        """
        Concentration fitting function from https://arxiv.org/pdf/1502.00391, this one assumes that the mass of the host halo is given in
        M200
        """
        if self.z <= 4:
            alpha = 1.7543 - 0.2766 * (1 + self.z) + 0.02039 * (1 + self.z) ** 2
            beta = 0.2753 + 0.00351 * (1 + self.z) - 0.3038 * (1 + self.z) ** 0.0269
            gamma = -0.01537 + 0.02102 * (1 + self.z) ** -0.1475
            return 10 ** (alpha + beta * self.logmass * (1 + gamma * (self.logmass) ** 2))
        elif self.z > 4:
            alpha = 1.3081 - 0.1078 * (1 + self.z) + 0.00398 * (1 + self.z) ** 2
            beta = 0.0223 - 0.0944 * (1 + self.z) ** -0.3907
            return 10 ** (alpha + beta * self.logmass)

    @staticmethod
    def f_NFW(c):
        return np.log1p(c) - c / (1 + c)


class Host_precomputed():
    """
    This is an alternative class to the standard Host_Halo class in Insert_Halo. Here we assume that the host halo lensing data has been computed/fit to 
    observational data. Doing this results in a lensmodel which we then use here to attempt to compute the host halo mass for the simulation.
    """
    def __init__(self, IHB, precomputed_host_data):
        self.IHB = IHB
        self._z = self.IHB.zlens
        self.precomputed_host_data = precomputed_host_data
        self._mass = self.estimate_mass_from_data()
        #self.profile = "NFW"
        #self.GM = geometry_and_mass_function
        #self.sh = self.GM.SH.sh_settings

    @property
    def logmass(self):
        return self._mass

    @property
    def mass(self):
        return 10 ** self._mass

    @property
    def z(self):
        return self._z


    def estimate_mass_from_data(self):
        """
        This will be difficult, but here the host halo mass should be computed from the provided lensing data (first find a paper on how to do this....)
        """
        pass















