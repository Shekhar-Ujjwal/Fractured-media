"""
HTI.py
------
Computes HTI stiffness coefficients from fracture weakness parameters.

Dependencies:
    numpy
    scipy.io

Author: Ujjwal Shekhar
"""

import numpy as np
from scipy.io import loadmat, savemat

#Case 1
def build_hti_model(fracture_model_mat="fracture_weakness_model.mat",
                    output_mat="HTI_model.mat"):
    """
    Build HTI stiffness coefficients from fracture weakness parameters.

    Parameters
    ----------
    fracture_model_mat : str
        Path to .mat file containing delN and delT.
    output_mat : str
        Output .mat file name for stiffness coefficients.

    Returns
    -------
    dict
        Dictionary containing HTI stiffness components.
    """
#Case 2
#def build_hti_model(fracture_model_mat="fracture_weakness_model.mat",
#    background_model_mat="background_stiffness.mat",
#    output_mat="HTI_model.mat"):
    
    # Load fracture weakness model
    data = loadmat(fracture_model_mat)
    delN = data["delN"]
    delT = data["delT"]

    # Physical parameters (homogeneous isotropic background)
    rho = 2500.0   # kg/m^3
    Vp  = 3231.0   # m/s
    Vs  = 1844.0   # m/s

    # Background stiffness coefficients
    c11_b = rho * Vp**2
    c33_b = c11_b
    c44_b = rho * Vs**2
    c66_b = c44_b
    c12_b = c11_b - 2 * c44_b
    c13_b = c12_b

    # HTI stiffness coefficients (vectorized)
    c11 = c11_b * (1 - delN)
    c12 = c12_b * (1 - delN)
    c13 = c13_b * (1 - delN)

    c22 = c11_b * (1 - delN * (c12_b / c11_b)**2)
    c23 = c13_b * (1 - delN * (c12_b / c11_b))

    c33 = c33_b * (1 - delN * ((c13_b**2) / (c33_b * c11_b)))

    c44 = c44_b * np.ones_like(delN)
    c55 = c44_b * (1 - delT)
    c66 = c44_b * (1 - delT)

    # Save to MATLAB-compatible file
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
    build_hti_model()
