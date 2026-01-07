# Elastic parameters for fractured-anisotropic media
Effective seismic models from fractured rocks

## Overview
This repository contains research code for computing **stiffness coefficients in fractured media**
from **fracture parameters and background elastic parameters**. The formulation is based on the **linear slip theory**.

The code is intended for **research and educational purposes** and is under active development.

---

## Features
- Useful for modelling seismic wave propagation in HTI, orthorhombic and monoclinic media
- Inputs fracture weaknesses (parameters accounting for fracture density and fracture-filling materials)
- Assumes thin vertical planar fractures
- Easy to implement in anisotropic-elastic full waveform inversion (FWI) framework

---

## Repository Structure
.
├── src/                # Source code
│   ├── main.py
├── README.md
├── LICENSE
├── requirements.txt
└── .gitignore

## References
1. Schoenberg, M. (1983) Reflection of elastic waves from periodically stratified media with interfacial slip, Geophysical Prospecting, 31(2), 265–292.
2. Schoenberg, M. and Helbig, K. (1997) Orthorhombic media: Modeling elastic wave behavior in a vertically fractured earth, Geophysics, 62(6), 1954–1974.
3. Shekhar, U. (2021) Effective seismic model from fractured rock. Master thesis. Trondheim, Norway: NTNU. Available at ntnuopen.ntnu.no/ntnu-xmlui/handle/11250/2834586.
4. Shekhar, U., Jakobsen, M., Pšenčík, I. and Xiang, K. (2025) Seismic full waveform inversion for fracture parameters in anisotropic media
   , Geophysical Prospecting, 73(5), 1606–1634, 2025.
