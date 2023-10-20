import numpy as np

class BuildingAnalysis:
    """
BuildingAnalysis Class:
-----------------------
This class represents a simplified analysis for a building subjected to wind loads. The analysis is based on predefined building dimensions, wind pressures, material properties of CLT (Cross Laminated Timber) panels, connectors, and a foundation. 

Attributes:
- Geometry: Includes attributes related to the building's core, overall dimensions, and the maximum height of a CLT panel.
- Forces: Contains attributes related to wind pressures and other floor forces.
- CLT material: Properties of a fictitious CLT panel built up of 7 layers each 40 mm thick.
- Connectors: Attributes related to shear connectors, hold down connectors, and pressure plates.
- Foundation: Includes a fictitious spring constant for the foundation.

Methods:
- compute_moments: Computes the moments at three different heights of the building due to wind loads.
- compute_shear_forces: Computes the shear forces at three different heights of the building due to wind loads.
- compute_displacements: Calculates various displacements in the building such as bending, shear, connector shear, hold down, and foundation displacements.
- check_displacements: Computes the total displacement of the building and checks it against a predefined limit based on building height.
- display_results: Prints out the computed moments, shear forces, individual displacements, total displacement, and displacement check values.

Usage:
The class can be instantiated and its methods called to conduct an analysis of the building under predefined conditions. Results can be displayed using the display_results method.

Note:
The class represents a simplified building analysis model, and actual real-world scenarios may require a more comprehensive analysis.
"""
    def __init__(self):
        # Geometry
        self.b_core = 6000  # mm 
        self.b_building = self.b_core*3  # mm
        self.h = 3500*12  # height of the building
        self.h_clt = 2*3500  # max height clt panel is 6m 

        # Forces
        self.p_wind = 1.25E-3  # Mpa
        self.q_wind = self.p_wind * self.b_building
        self.G_floor = self.b_core * self.b_core * 2.0E-3
        self.Q_floor = self.b_core * self.b_core * 2.5E-3
        self.Q_floor_c = self.Q_floor * 0.4

        # CLT material, ficticious buildup 7 layers 40 mm
        self.t_0 = 4*40
        self.G_m = 690  # N/mm2
        self.G = 0.75*self.G_m
        self.E_0 = 11000  # N/mm2
        self.I_net = 2*((1/12) * self.t_0*self.b_core**3 + (1/12)*self.b_core*self.t_0**3 + (self.t_0*self.b_core)*(self.b_core / 2)**2)
        self.A_s = 2*self.b_core*self.t_0

        # Connectors
        self.K_shear = 8.2E3
        self.K_holddown = 13E3
        self.K_pplate = self.K_holddown
        self.n_shear = 8
        self.n_holddown = 6

        # foundation
        self.k_r = 10E100

    def compute_moments(self):
        M_0 = self.q_wind * (self.h**2) * (1/2)
        M_7 = self.q_wind * ((self.h*(2/3))**2) * (1/2)
        M_14 = self.q_wind * ((self.h*(1/3))**2) * (1/2)
        return [M_0, M_7, M_14]

    def compute_shear_forces(self):
        V_0 = self.q_wind * self.h
        V_7 = self.q_wind * self.h * (2/3)
        V_14 = self.q_wind * self.h * (1/3)
        return [V_0, V_7, V_14]

    def compute_displacements(self):
        u_b = self.q_wind * ((self.h)**4) / (8 * self.E_0 * self.I_net)
        u_s = (self.q_wind * self.h**2) / (2 * self.G * self.A_s)
        u_f = (self.q_wind * self.h**3) / (2 * self.k_r)

        F_cs = [x / (2 * self.n_shear) for x in self.compute_shear_forces()]
        u_cs = [x / self.K_shear for x in F_cs]
        u_cs_tot = np.sum(u_cs)

        F_hd = [x / (self.b_core * self.n_holddown) for x in self.compute_moments()]
        d_hd = [x / self.K_holddown for x in F_hd]
        theta = [x / self.b_core for x in d_hd]
        u_ch = [x * self.h_clt for x in theta]
        u_ch_tot = np.sum(u_ch)

        return [u_b, u_s, u_cs_tot, u_ch_tot, u_f]

    def check_displacements(self):
        u_list = self.compute_displacements()
        u_tot = np.sum(u_list)
        u_check = self.h / 500
        return u_tot, u_check

    def display_results(self):
        M = self.compute_moments()
        V = self.compute_shear_forces()
        u_list = self.compute_displacements()
        u_tot, u_check = self.check_displacements()

        print("Moments:", M)
        print("Shear Forces:", V)
        print("Displacement List:", u_list)
        print("Total Displacement:", u_tot)
        print("Displacement Check:", u_check)

if __name__ == '__main__':
    analysis = BuildingAnalysis()
    analysis.display_results()
