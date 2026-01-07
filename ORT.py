"""
VFTI.py
-------
Code to generate stiffness coefficients for heterogeneous VFTI /
orthorhombic media based on linear slip theory (Schoenberg, 1983)

Author: Ujjwal Shekhar 
(ChatGPT has been used for polishing the code, notations are self-explanatory, ask author if not clear or check MATLAB code)
"""

import numpy as np
from scipy.io import loadmat, savemat


def build_vfti_model(
    fracture_model_mat="fracture_weakness_model.mat",
    density_model_mat=None,
    velocity_model_mat=None,
    output_mat="VFTI_model.mat"
):
    """
    Builds VFTI / orthorhombic stiffness coefficients.

    If density_model_mat and velocity_model_mat are None:
        -> Case 1: Homogeneous VTI background
    Else:
        -> Case 2: Heterogeneous VTI background
    """

    # ------------------------------------------------------------
    # Load fracture weaknesses
    frac = loadmat(fracture_model_mat)
    delN = frac["delN"]
    delV = frac["delV"]
    delH = frac["delH"]

    # ------------------------------------------------------------
    # CASE 1: Homogeneous VTI background
    if density_model_mat is None and velocity_model_mat is None:

        rho = 2500  # kg/m^3
        Vp = 3231   # m/s
        Vs = 1844.  # m/s

        c11_0 = rho * Vp**2
        c44_0 = rho * Vs**2
        c12_0 = c11_0 - 2 * c44_0

        C0 = np.array([
            [c11_0, c12_0, c12_0, 0,      0,      0],
            [c12_0, c11_0, c12_0, 0,      0,      0],
            [c12_0, c12_0, c11_0, 0,      0,      0],
            [0,      0,      0,     c44_0, 0,      0],
            [0,      0,      0,     0,      c44_0, 0],
            [0,      0,      0,     0,      0,      c44_0]
        ])

        # A-parameters (Pšenčík & Farra, 2024)
        eps_x = 0.103
        eps_z = 0.0
        gamma_x = 0.0
        gamma_z = 0.106
        eta = -0.283

        c11_b = C0[0, 0] * (1.0 + 2.0 * eps_x)
        c12_b = C0[0, 0] * (1.0 + 2.0 * eps_x) - 2.0 * C0[3, 3] * (1.0 + 2.0 * gamma_z)
        c13_b = C0[0, 0] * (1.0 + eps_x + eps_z + eta) - 2.0 * C0[3, 3] * (1.0 + 2.0 * gamma_x)
        c33_b = C0[0, 0] * (1.0 + 2.0 * eps_z)
        c44_b = C0[3, 3] * (1.0 + 2.0 * gamma_x)
        c66_b = C0[3, 3] * (1.0 + 2.0 * gamma_z)

    # ------------------------------------------------------------
    # CASE 2: Heterogeneous VTI background
    else:
        rho = loadmat(density_model_mat)["rho"]

        vel = loadmat(velocity_model_mat)
        c11_b = vel["c11_b"]
        c33_b = vel["c33_b"]
        c44_b = vel["c44_b"]
        c66_b = vel["c66_b"]
        c13_b = vel["c13_b"]

        c12_b = c11_b - 2.0 * c66_b

    # ------------------------------------------------------------
    # Elastic parameters of VFTI / orthorhombic media
    c11 = c11_b * (1 - delN)
    c12 = c12_b * (1 - delN)
    c13 = c13_b * (1 - delN)

    c22 = c11_b * (1 - delN * (c12_b / c11_b)**2)
    c23 = c13_b * (1 - delN * (c12_b / c11_b))

    c33 = c33_b * (1 - delN * ((c13_b**2) / (c33_b * c11_b)))

    c44 = c44_b * np.ones_like(delN)
    c55 = c44_b * (1 - delV)
    c66 = c66_b * (1 - delH)

    # ------------------------------------------------------------
    # Save model
    savemat(
        output_mat,
        {
            "c11": c11,
            "c12": c12,
            "c13": c13,
            "c22": c22,
            "c23": c23,
            "c33": c33,
            "c44": c44,
            "c55": c55,
            "c66": c66,
        }
    )

    return {
        "c11": c11,
        "c12": c12,
        "c13": c13,
        "c22": c22,
        "c23": c23,
        "c33": c33,
        "c44": c44,
        "c55": c55,
        "c66": c66,
    }


if __name__ == "__main__":
    build_vfti_model()
