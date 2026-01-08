"""
Effective monoclinic model for two non-orthogonal vertical fracture sets
in a VTI background

Author: Ujjwal Shekhar 
(ChatGPT has been used for polishing the code, notations are self-explanatory, ask author if not clear)
"""

import numpy as np


def monoclinic_two_fracture_sets():
    # ------------------------------------------------------------
    # Fracture weaknesses (normal, vertical-tangential, horizontal-tangential)
    delN1, delV1, delH1 = 0.10, 0.15, 0.20
    delN2, delV2, delH2 = 0.06, 0.10, 0.15

    # Fracture azimuths (radians)
    phi1 = np.pi / 4.0
    phi2 = np.pi / 2.0

    # ------------------------------------------------------------
    # VTI background parameters
    rho = 2300.0          # kg/m^3
    Vp0 = 2600.0          # m/s (P-wave vertical velocity)
    Vs0 = 1200.0          # m/s (S-wave vertical velocity)

    # Background stiffnesses (GPa)
    C33b = rho * Vp0**2 / 1e9
    C44b = rho * Vs0**2 / 1e9
    
    # Thomsen (1986) parameters
    epsilon = 0.1
    gamma = 0.15
    delta = 0.05

    C11b = (2.0 * epsilon + 1.0) * C33b
    C66b = (2.0 * gamma + 1.0) * C44b
    C13b = np.sqrt(
        2 * C33b * (C33b - C44b) * delta + (C33b - C44b)**2
    ) - C44b
    C12b = C11b - 2.0 * C66b

    # ------------------------------------------------------------
    # VTI background stiffness matrix
    Cb = np.array([
        [C11b, C12b, C13b, 0,     0,     0],
        [C12b, C11b, C13b, 0,     0,     0],
        [C13b, C13b, C33b, 0,     0,     0],
        [0,     0,     0,     C44b, 0,     0],
        [0,     0,     0,     0,     C44b, 0],
        [0,     0,     0,     0,     0,     C66b]
    ])

    # ------------------------------------------------------------
    # Fracture compliance matrices
    Kn1 = delN1 / (C11b * (1.0 - delN1))
    Kv1 = delV1 / (C44b * (1.0 - delV1))
    Kh1 = delH1 / (C66b * (1.0 - delH1))

    Kn2 = delN2 / (C11b * (1.0 - delN2))
    Kv2 = delV2 / (C44b * (1.0 - delV2))
    Kh2 = delH2 / (C66b * (1.0 - delH2))

    D1 = np.array([
        [Kn1, 0,   0,   0,   0,    0],
        [0,   0,   0,   0,   0,    0],
        [0,   0,   0,   0,   0,    0],
        [0,   0,   0,   0,   0,    0],
        [0,   0,   0,   0,   Kv1,  0],
        [0,   0,   0,   0,   0,    Kh1]
    ])

    D2 = np.array([
        [Kn2, 0,   0,   0,   0,    0],
        [0,   0,   0,   0,   0,    0],
        [0,   0,   0,   0,   0,    0],
        [0,   0,   0,   0,   0,    0],
        [0,   0,   0,   0,   Kv2,  0],
        [0,   0,   0,   0,   0,    Kh2]
    ])

    # ------------------------------------------------------------
    # Bond rotation matrix
    def rotation_matrix(phi):
        return np.array([
            [np.cos(phi)**2,         np.sin(phi)**2,        0, 0, 0,  np.sin(2*phi)],
            [np.sin(phi)**2,         np.cos(phi)**2,        0, 0, 0, -np.sin(2*phi)],
            [0,                       0,                    1, 0, 0,  0],
            [0,                       0,                    0,  np.cos(phi), -np.sin(phi), 0],
            [0,                       0,                    0,  np.sin(phi),  np.cos(phi), 0],
            [-0.5*np.sin(2*phi),      0.5*np.sin(2*phi),    0, 0, 0,  np.cos(2*phi)]
        ])

    R1 = rotation_matrix(phi1)
    R2 = rotation_matrix(phi2)

    # ------------------------------------------------------------
    # Monoclinic stiffness matrix (GPa)
    I6 = np.eye(6)
    Cmono = Cb @ np.linalg.inv(
        I6 + (R1 @ D1 @ R1.T + R2 @ D2 @ R2.T) @ Cb
    )
    # ------------------------------------------------------------
    # Density-normalized stiffness coefficients
    A = {
        "A11": Cmono[0, 0] / rho,
        "A12": Cmono[0, 1] / rho,
        "A13": Cmono[0, 2] / rho,
        "A16": Cmono[0, 5] / rho,
        "A22": Cmono[1, 1] / rho,
        "A23": Cmono[1, 2] / rho,
        "A26": Cmono[1, 5] / rho,
        "A33": Cmono[2, 2] / rho,
        "A36": Cmono[2, 5] / rho,
        "A44": Cmono[3, 3] / rho,
        "A45": Cmono[3, 4] / rho,
        "A55": Cmono[4, 4] / rho,
        "A66": Cmono[5, 5] / rho,
    }

    return Cmono, A
    
# Note that any number of fracture sets can be included 
# by adding their compliances in a similar fashion

if __name__ == "__main__":
    Cmono, A = monoclinic_two_fracture_sets()
    print("Monoclinic stiffness matrix (GPa):")
    print(Cmono)
    print("\\nDensity-normalized stiffness coefficients:")
    for k, v in A.items():
        print(f"{k} = {v:.6e}")
