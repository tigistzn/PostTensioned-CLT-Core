from scipy.optimize import minimize
import pandas as pd

class PretensionedCLT:
    """
    This class facilitates the analysis of pretensioned Cross Laminated Timber (CLT) plates.
    
    Attributes:
    - Various physical and mechanical properties of the CLT, such as areas, moduli of elasticity,
      stresses, and geometric properties.
    - Limit state values including moments, forces, and displacements.
    
    Methods:
    - Initial computations for the CLT properties.
    - Calculate base values related to the axial force due to self-weight.
    - Assess limit states like decompression, elastic, and yield, among others.
    - Determine the compression zone width based on specific states.
    - Calculate bending displacements for various states.
    - Export the calculated results to an Excel file.
    """
    
    def __init__(self, A_pt=1050, E_pt=205E3, f_pi=558, gam=5.5, h=4000, lw=1400, f_c0=24, t_w=150, tl=10, EI_ef=10E3 * 22275E4):
        """
        Initialize the CLT properties. Default values are used if no arguments are provided.

        :param A_pt: Int or Float, cross-sectional area, default is 1050.
        :param E_pt: Int or Float, Elastic modulus, default is 205E3.
        :param f_pi: Int or Float, stress, default is 558.
        :param gam: Int or Float, weight per unit volume, default is 5.5.
        :param h: Int, height, default is 4000.
        :param lw: Int, width, default is 1400.
        :param f_c0: Int or Float, "yield" stress of CLT, default is 24.
        :param t_w: Int, thickness of the CLT, default is 150.
        :param tl: Int, tolerance for compression zone check, default is 10.
        :param EI_ef: Int or Float, effective flexural stiffness, default is 10E3 * 22275E4.
        """
        self.A_pt = A_pt
        self.E_pt = E_pt
        self.f_pi = f_pi
        self.F_pi = A_pt * f_pi
        self.gam = gam
        self.h = h
        self.lw = lw
        self.f_c0 = f_c0
        self.t_w = t_w
        self.H_cr = 2 * t_w
        self.H_pt = h
        self.eps_c0 = 0.008
        self.eps_cs = 0.02
        self.eps_cu = 0.05
        self.tl = tl
        self.EI_ef = EI_ef
        self.eps_pi = self.F_pi / E_pt
        
    def calculate_base_values(self):
        """Calculate the basic values like axial force due to self-weight."""
        self.N_g = (self.h/1000 * self.lw/1000) * (self.t_w/1000) * self.gam * 1000

    def calculate_limit_states(self):
        """Calculate actions (moment and shear) for different limit states."""
        # DEC - Decompression limit state
        self.C_dec = self.N_g + self.F_pi
        self.M_dec = ((self.F_pi + self.N_g) * (2/3) - (self.F_pi + self.N_g) / 2) * self.lw
        self.V_dec = self.M_dec / self.h

        # ELL - Elastic limit state
        self.M_ell = (self.N_g + self.F_pi) * (1/2 - 3/8*1/3) * self.lw
        self.V_ell = self.M_ell / self.h

        # YCLT - yielding of the CLT
        self.c_yclt = self._compression_zone_width(self.eps_c0, self.CH_yclt)
        C_yclt = 1/2 * self.f_c0 * self.c_yclt * self.t_w
        self.M_yclt = (C_yclt * (self.lw - (1/3) * self.c_yclt)) - (self.N_g + self.F_pi) * (1/2) * self.lw
        self.V_yclt = self.M_yclt / self.h

        # SCLT
        self.c_sclt = self._compression_zone_width(self.eps_cs, self.CH_sclt)
        C_sclt = ((1 - self.eps_c0 / self.eps_cs) + (1/2) * (self.eps_c0 / self.eps_cs)) * self.c_sclt * self.f_c0 * self.t_w
        self.M_sclt = (C_sclt * (self.lw - (1/3) * self.c_sclt)) - (self.N_g + self.F_pi) * (1/2) * self.lw
        self.V_sclt = self.M_sclt / self.h

        # CCLT
        self.c_cclt = self._compression_zone_width(self.eps_cu, self.CH_cclt)
        C_cclt = ((1 - self.eps_c0 / self.eps_cu) + (1/2) * (self.eps_c0 / self.eps_cu)) * self.c_cclt * self.f_c0 * self.t_w
        self.M_cclt = (C_cclt * (self.lw - (1/3) * self.c_cclt)) - (self.N_g + self.F_pi) * (1/2) * self.lw
        self.V_cclt = self.M_cclt / self.h
        
    def _compression_zone_width(self, eps, ch_func):
        """Find the compression zone width using the given check function."""
        return minimize(ch_func, 300, tol=self.tl).x[0]
    
    def CH_yclt(self, c_yclt):
        """
        Compute the check function for the yielding of the CLT.
        
        Parameters:
        - c_yclt (float): Compression zone width for YCLT.
        
        Returns:
        - float: Value of the check function.
        """
        C_yclt = 0.5 * self.f_c0 * c_yclt * self.t_w 
        eps_p_yclt= self.eps_pi + (self.H_cr / self.H_pt) * ((0.5 * self.lw - c_yclt) / c_yclt) * self.eps_yclt
        F_p_yclt = eps_p_yclt * self.E_pt
        Check_yclt = abs(F_p_yclt + self.N_g - C_yclt)
        return Check_yclt

    def CH_sclt(self, c_sclt):
        """
        Compute the check function for the SCLT state.
        
        Parameters:
        - c_sclt (float): Compression zone width for SCLT.
        
        Returns:
        - float: Value of the check function.
        """
        C_sclt = ((1 - self.eps_c0 / self.eps_cs) + 0.5 * (self.eps_c0 / self.eps_cs)) * c_sclt * self.f_c0 * self.t_w
        eps_p_sclt = self.eps_pi + (self.H_cr / self.H_pt) * ((0.5 * self.lw - c_sclt) / c_sclt) * self.eps_sclt
        F_p_sclt = eps_p_sclt * self.E_pt
        Check_sclt = abs(F_p_sclt + self.N_g - C_sclt)
        return Check_sclt

    def CH_cclt(self, c_cclt):
        """
        Compute the check function for the CCLT state.
        
        Parameters:
        - c_cclt (float): Compression zone width for CCLT.
        
        Returns:
        - float: Value of the check function.
        """
        eps_p_cclt= self.eps_pi + (self.H_cr / self.H_pt) * ((0.5 * self.lw - c_cclt) / c_cclt) * self.eps_cclt
        C_cclt = ((1 - self.eps_c0 / self.eps_cclt) + 0.5 * (self.eps_c0 / self.eps_cclt)) * c_cclt * self.f_c0 * self.t_w
        F_p_cclt = eps_p_cclt * self.E_pt
        Check_cclt = abs(F_p_cclt + self.N_g - C_cclt)
        return Check_cclt

    def displacements(self):
        """Calculate bending displacements for each limit state."""
        self.ub_dec = self._ubending(self.V_dec)
        self.ub_yclt = self._ubending(self.V_yclt)
        self.ub_sclt = self._ubending(self.V_sclt)
        self.ub_cclt = self._ubending(self.V_cclt)
        
    def _ubending(self, F):
        """Compute the bending displacement."""
        u = (F * self.h**3) / (3 * self.EI_ef)
        return u
    
    def write_to_excel(self, filename="results.xlsx"):
        """Save the calculated values to an Excel file."""
        
        # Create a dataframe to hold the results
        df = pd.DataFrame({
            'Parameter': ['Base axial force due to self-weight', 
                          'DEC Compressive force', 'DEC Maximum moment', 'DEC Shear force', 'DEC Displacement',
                          'ELL Maximum moment', 'ELL Shear force',
                          'YCLT Compression width', 'YCLT Maximum moment', 'YCLT Shear force', 'YCLT Displacement',
                          'SCLT Compression width', 'SCLT Maximum moment', 'SCLT Shear force', 'SCLT Displacement',
                          'CCLT Compression width', 'CCLT Maximum moment', 'CCLT Shear force', 'CCLT Displacement'],
            'Value': [self.N_g, 
                      self.C_dec, self.M_dec, self.V_dec, self.ub_dec,
                      self.M_ell, self.V_ell,
                      self.c_yclt, self.M_yclt, self.V_yclt, self.ub_yclt,
                      self.c_sclt, self.M_sclt, self.V_sclt, self.ub_sclt,
                      self.c_cclt, self.M_cclt, self.V_cclt, self.ub_cclt],
        })

        # Save the dataframe to Excel
        with pd.ExcelWriter(filename, engine='openpyxl') as writer:
            df.to_excel(writer, index=False)

# Using the class
clt = PretensionedCLT()
clt.calculate_base_values()
clt.calculate_limit_states()
clt.displacements()
clt.write_to_excel()