import math

class BuildingGeometry:
    """
    Represents the geometrical properties of a building.
    
    Attributes:
        b_core (float): Core breadth of the building.
        b_building (float): Total breadth of the building.
        h (float): Total height of the building.
        A_floor (float): Area of a single floor.
    """
    def __init__(self, b_core, num_floors, floor_height, floor_length):
        self.b_core = b_core
        self.b_building = self.b_core * 3
        self.h = floor_height * num_floors  # height of the building
        self.A_floor = floor_length ** 2  # Area of the building floor

class CLTMaterial:
    """
    Represents the material properties of the Cross Laminated Timber (CLT) used in the building.
    
    Attributes:
        geometry (BuildingGeometry): The geometry of the building.
        t_0 (float): Thickness of the material.
        G_m (float): Shear modulus of the material.
        G (float): Adjusted shear modulus.
        E_0 (float): Elastic modulus of the material.
        W (float): Weight of the material.
        f_c0 (float): Compression property of the material.
    """
    def __init__(self, building_geometry, t_0, G_m, E_0, W, f_c0):
        self.geometry = building_geometry
        self.t_0 = t_0
        self.G_m = G_m
        self.G = 0.75 * self.G_m
        self.E_0 = E_0
        self.W = W
        self.f_c0 = f_c0

class Forces:
    """
    Represents various forces that act on the building.
    
    Attributes:
        geometry (BuildingGeometry): The geometry of the building.
        p_wind (float): Wind pressure on the building.
        q_wind (float): Wind force on the building.
        G_floor (float): Dead load on a floor.
        Q_floor (float): Live load on a floor.
        psi (float): Load factor.
        f (float): Load duration factor.
    """
    def __init__(self, building_geometry, p_wind, G_floor, Q_floor, psi, f):
        self.geometry = building_geometry
        self.p_wind = p_wind
        self.q_wind = self.p_wind * self.geometry.b_building
        self.G_floor = G_floor
        self.Q_floor = Q_floor
        self.psi = psi
        self.f = f

class Connectors:
    """
    Represents various connectors used in post-tensioning.
    
    Attributes:
        K_shear (float): Shear stiffness.
        K_holddown (float): Holddown stiffness.
        K_pplate (float): Plate stiffness.
        n_shear (int): Number of shear connectors.
        n_holddown (int): Number of holddown connectors.
    """
    def __init__(self, K_shear, K_holddown, K_pplate, n_shear, n_holddown):
        self.K_shear = K_shear
        self.K_holddown = K_holddown
        self.K_pplate = K_pplate
        self.n_shear = n_shear
        self.n_holddown = n_holddown

class PostTensionLoss:
    """
    Calculates and represents post-tension losses in the building system.
    
    Attributes:
        material (CLTMaterial): Material properties.
        k (float): Efficiency factor.
        fpk (float): Prestress force.
        A_section (float): Cross-sectional area.
        forces (Forces): Forces acting on the building.
        connectors (Connectors): Connector properties.
    """
    def __init__(self, material, k, fpk, A_section, forces, connectors):
        self.material = material
        self.k = k
        self.fpk = fpk
        self.A_section = A_section
        self.forces = forces
        self.connectors = connectors

        self.compute_values()

    def compute_values(self):
        """Calculates various post-tension related values."""
        self.N_g_k = self.material.W * self.material.geometry.A_floor * self.material.geometry.h
        self.N_min_k = 6 * self.material.geometry.A_floor * self.forces.G_floor
        self.N_max_k = 6 * self.material.geometry.A_floor * (self.forces.G_floor + self.forces.Q_floor)

        self.N_G_d = self.N_g_k * 1.2
        self.N_min_d = self.N_min_k * 1.2 * 0.9
        self.N_max_d = 1.2 * 6 * self.material.geometry.A_floor * self.forces.G_floor + 1.5 * 6 * self.material.geometry.A_floor * self.forces.Q_floor

        self.M_wind_k = (1/2) * self.forces.q_wind * self.material.geometry.h**2
        self.M_wind_d = 1.5 * self.M_wind_k

        I_net = 2 * ((1/12) * self.material.t_0 * self.material.geometry.b_core**3 + (1/12) * self.material.geometry.b_core * self.material.t_0**3 + (self.material.t_0 * self.material.geometry.b_core) * (self.material.geometry.b_core / 2)**2)
        A_net = 4 * self.material.geometry.b_core * self.material.t_0

        sigma_m_d = (self.M_wind_d * (self.material.geometry.b_building / 2)) / I_net
        sigma_nmin_d = self.N_min_d / A_net
        sigma_nmax_d = self.N_max_d / A_net

        P_min = A_net * (sigma_m_d - sigma_nmin_d)
        P_max = A_net * (self.material.f_c0 - sigma_m_d - sigma_nmax_d)

        self.P0 = P_min / self.k
        A_bar = (1/4) * math.pi * 50**2
        self.P_bar = 0.3 * self.fpk * A_bar

        # Vertical shortening due to the prestress
        self.d_v = (self.P0 * self.material.geometry.h) / (self.material.E_0 * A_net)

    def get_results(self):
        """Returns vertical shortening and prestress loss values."""

        return self.d_v, self.P0

def main():
    """
    Entry point of the program to compute and display the post-tension losses for a CLT building.

    Initializes the building geometry, material properties, external forces, and connectors. 
    Calculates the post-tension losses based on these inputs and then displays 
    the vertical shortening and prestress loss values.
    """
    geometry = BuildingGeometry(6000, 6, 3500, 12000)
    material = CLTMaterial(geometry, 160, 690, 11000, 5.5E-9, 13.4)
    forces = Forces(geometry, 1.25E-3, 2.0E-3, 2.5E-3, 0.4, 0.9)
    connectors = Connectors(8.2E3, 13E3, 13E3, 8, 6)
    post_tension = PostTensionLoss(material, 0.5, 1230, 4*280*6000, forces, connectors)

    vertical_shortening, prestress_loss = post_tension.get_results()
    print(f"Vertical Shortening: {vertical_shortening} mm")
    print(f"Prestress Loss: {prestress_loss} N")

if __name__ == '__main__':
    main()