import numpy as np
import pandas as pd
from scipy import optimize

class Wall:
    """
    A class to represent the CLT wall and the forces acting on the wall. 
    """
    def __init__(self, total_height, wall_width, rod_eccentricity, wall_thickness, normal_force, prestress_force):
        """
        Initialize the Wall instance with given parameters.

        Parameters:
        - total_height (float): Total height of the wall.
        - wall_width (float): Width of the wall.
        - rod_eccentricity (float): Eccentricity of the rod.
        - wall_thickness (float): Thickness of the wall.
        - normal_force (float): The normal force applied to the wall.
        - prestress_force (float): The prestressing force applied to the wall.
        """
        self.TOTAL_HEIGHT = total_height
        self.WALL_WIDTH = wall_width
        self.ROD_ECCENTRICITY = rod_eccentricity
        self.WALL_THICKNESS = wall_thickness
        self.ROD_HEIGHT = total_height
        self.NORMAL_FORCE = normal_force
        self.YOUNGS_MODULULS = 10E3
        self.INERTIA = 3 * (1/12) * 30 * self.WALL_WIDTH**3
        self.GA_VALUE = 0.75 * 690 * self.WALL_WIDTH * self.WALL_THICKNESS
        self.EPS_YCLT = 0.008
        self.EPS_CCLT = 0.02
        self.PRESTRESS_FORCE = prestress_force
        self.E_P = 205E3
        self.AREA_P = 0.25 * np.pi * 24**2

    def bending(self, M):
        """calculate the bending deformation"""
        q = 2 * M / (self.TOTAL_HEIGHT**2)
        u = q * self.TOTAL_HEIGHT**4 / (8 * self.YOUNGS_MODULULS * self.INERTIA)
        return u

    def shear(self, M):
        """calculate the shear deformation"""
        q = 2 * M / (self.TOTAL_HEIGHT**2)
        u = q * self.TOTAL_HEIGHT**2 / (2 * self.GA_VALUE)
        return u

    def FindHeight_1(self, c):
        """
        Calculate the difference between applied forces and the normal force for the yielding case.

        Parameter:
        - c (float): Distance from the neutral axis.

        Returns:
        - float: Difference between forces.
        """
        H_cr = 2 * self.WALL_THICKNESS
        l_pt = self.WALL_WIDTH - self.ROD_ECCENTRICITY
        deps_p = (H_cr / self.ROD_HEIGHT) * ((0.5 * l_pt - c) / c) * self.EPS_YCLT
        eps_p = (self.PRESTRESS_FORCE / 2) / (self.AREA_P * self.E_P) + deps_p
        F_p = eps_p * self.E_P * self.AREA_P
        F_clt = 0.5 * 24 * c * self.WALL_THICKNESS
        return abs(F_p - F_clt - self.NORMAL_FORCE)

    def FindHeight_2(self, c):
        """
        Calculate the difference between applied forces and the normal force for the crushing case.

        Parameter:
        - c (float): Distance from the neutral axis.

        Returns:
        - float: Difference between forces.
        """
        H_cr = 2 * self.WALL_THICKNESS
        l_pt = self.WALL_WIDTH - self.ROD_ECCENTRICITY
        deps_p = (H_cr / self.ROD_HEIGHT) * ((0.5 * l_pt - c) / c) * self.EPS_CCLT
        eps_p = (self.PRESTRESS_FORCE / 2) / (self.AREA_P * self.E_P) + deps_p
        F_p = eps_p * self.E_P * self.AREA_P
        F_c = 0.5 * self.EPS_YCLT / self.EPS_CCLT * c * self.WALL_THICKNESS * 24 + (1 - self.EPS_YCLT / self.EPS_CCLT) * c * self.WALL_THICKNESS * 24
        return abs(F_p - F_c - self.NORMAL_FORCE)
    
    def results(self):
        """
        Calculate results for various limits like decompression, elastic, yielding, and crushing.

        Returns:
        - dict: Dictionary containing displacement and stress values for each limit.
        """
        # Decompression limit calculations
        F_dec = self.F_pi + self.N
        M_dec = (1/6) * self.b_wall * F_dec
        u_dec_EI = self.bending(M_dec)
        u_dec_GA = self.shear(M_dec)
        u_dec_tot = u_dec_EI + u_dec_GA
        sigma_dec_p= (self.F_pi/2)/self.A_p
        sigma_dec_clt= F_dec/(self.t_wall * self.b_wall)

        # Elastic limit calculations
        F_ell = self.F_pi + self.N
        M_ell = (1/2-(3/8 * 1/3))*self.b_wall * F_ell
        u_ell_EI = self.bending(M_ell)
        u_ell_GA = self.shear(M_ell)
        u_ell_tot = u_ell_EI + u_ell_GA
        sigma_ell_p = (self.F_pi/2)/self.A_p
        sigma_ell_clt = F_ell/(self.t_wall*(3/8)*self.b_wall)

        # CLT yielding limit calculations
        c_yclt = optimize.minimize(self.FindHeight_1, 200, tol=1).x[0]
        theta = self.H_cr*self.eps_yclt/c_yclt 
        F_clt = 0.5*self.f_c0*c_yclt*self.t_eff
        diff_eps_p = (self.H_cr/self.H_pt) * ((0.5*self.b_wall-c_yclt)/c_yclt)*self.eps_yclt
        d_F_p = diff_eps_p*self.E_p*self.A_p
        M_yclt = d_F_p*(0.5*self.b_wall - self.d_s) + F_clt*(0.5*self.b_wall - (1/3)*c_yclt)
        u_yclt_EI = self.bending(M_yclt)
        u_yclt_GA = self.shear(M_yclt)
        u_yclt_theta = theta*self.h_total
        u_yclt_tot = u_yclt_EI + u_yclt_GA + u_yclt_theta

        sigma_yclt_p = self.F_pi/self.A_p + diff_eps_p*self.E_p
        sigma_yclt_clt = self.f_c0

        # CLT crushing limit calculations
        c_cclt = optimize.minimize(self.FindHeight_2, 200, tol=1).x[0]
        theta_cclt = self.H_cr*self.eps_cclt/c_yclt 
        F_cclt = 0.5*(self.eps_yclt/self.eps_cclt)*c_cclt*self.t_eff*self.f_c0 + (1-self.eps_yclt/self.eps_cclt)*c_cclt*self.t_eff*self.f_c0
        diff_eps_p = (self.H_cr/self.H_pt) * ((0.5*self.b_wall-c_cclt)/c_cclt)*self.eps_cclt
        d_F_p_cclt = diff_eps_p*self.E_p*self.A_p
        M_cclt = d_F_p_cclt*(0.5*self.b_wall - self.d_s) + F_cclt*(0.5*self.b_wall - (1/3)*c_cclt)
        u_cclt_EI = self.bending(M_cclt)
        u_cclt_GA = self.shear(M_cclt)
        u_cclt_theta = theta_cclt*self.h_total
        u_cclt_tot = u_cclt_EI + u_cclt_GA + u_cclt_theta
        sigma_cclt_p = self.F_pi/self.A_p + diff_eps_p*self.E_p
        sigma_cclt_clt = self.f_c0

        return {
            "Decompression Limit Displacement": u_dec_tot,
            "Elastic Limit Displacement": u_ell_tot,
            "CLT Yielding Limit Displacement": u_yclt_tot,
            "CLT Crushing Limit Displacement": u_cclt_tot,
            "Decompression Limit Stress (Steel)": sigma_dec_p,
            "Elastic Limit Stress (Steel)": sigma_ell_p,
            "CLT Yielding Limit Stress (Steel)": sigma_yclt_p,
            "CLT Crushing Limit Stress (Steel)": sigma_cclt_p,
            "Decompression Limit Stress (CLT)": sigma_dec_clt,
            "Elastic Limit Stress (CLT)": sigma_ell_clt,
            "CLT Yielding Limit Stress (CLT)": sigma_yclt_clt,
            "CLT Crushing Limit Stress (CLT)": sigma_cclt_clt,
        }
    
    def write_to_excel(self, filename="results.xlsx"):
        """Save the calculated values to an Excel file."""
        # Create a dataframe to hold the results
        df = pd.DataFrame({
            'Parameter': [
                'Decompression Limit Displacement', 'Elastic Limit Displacement', 'CLT Yielding Limit Displacement',
                'CLT Crushing Limit Displacement', 'Decompression Limit Stress (Steel)', 'Elastic Limit Stress (Steel)',
                'CLT Yielding Limit Stress (Steel)', 'CLT Crushing Limit Stress (Steel)', 'Decompression Limit Stress (CLT)',
                'Elastic Limit Stress (CLT)', 'CLT Yielding Limit Stress (CLT)', 'CLT Crushing Limit Stress (CLT)'
            ],
            'Value': [
                self.results()["Decompression Limit Displacement"], self.results()["Elastic Limit Displacement"],
                self.results()["CLT Yielding Limit Displacement"], self.results()["CLT Crushing Limit Displacement"],
                self.results()["Decompression Limit Stress (Steel)"], self.results()["Elastic Limit Stress (Steel)"],
                self.results()["CLT Yielding Limit Stress (Steel)"], self.results()["CLT Crushing Limit Stress (Steel)"],
                self.results()["Decompression Limit Stress (CLT)"], self.results()["Elastic Limit Stress (CLT)"],
                self.results()["CLT Yielding Limit Stress (CLT)"], self.results()["CLT Crushing Limit Stress (CLT)"]
            ],
        })
        with pd.ExcelWriter(filename, engine='openpyxl') as writer:
            df.to_excel(writer, index=False)
    
def main():
    """
    Main function for the Limnologen project case study to compute limit states for a CLT wall.

    In the context of the Limnologen project, this function simulates the behavior of a CLT wall subjected to 
    different loading conditions. The function initializes the Wall class based on given wall properties and loads. 
    It then computes various limit states, including decompression, elastic, CLT yielding, and CLT crushing. 
    The results of these calculations are saved to an Excel file named "results.xlsx".

    Parameters:
    None

    Returns:
    None


    Notes:
    - The function assumes specific wall dimensions, materials, and loadings as relevant to the Limnologen project.
    """
    total_height = 24500
    wall_width = 3700
    rod_eccentricity = 200
    wall_thickness = 5*30
    normal_force = 0
    prestress_force = 165E3*2

    wall_instance = Wall(total_height, wall_width, rod_eccentricity, wall_thickness, normal_force, prestress_force)
    wall_instance.write_to_excel("limnologen_wall_results.xlsx")
    
if __name__ == "__main__":
    main()
