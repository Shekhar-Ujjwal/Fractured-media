%  Code to generate stiffness coefficients for the heterogeneous VFTI (vertically fractured transversely isotropic media)
%  /orthorhombic model based on the linear slip theory (Schoenberg, 1983)
% Author: Ujjwal Shekhar

% Input the three fracture weaknesses at grid points 
load fracture_weakness_model delN delV delH

% delN = normal fracture weakness, delV = vertical-tangential fracture weakness, delH = horizontal-tangential fracture weakness

[Nx, Ny, Nz] = size(delN);  % Extracting the size of the model
 
%% Case 1: Homogeneous VTI background (in which vertical fractures are embedded):
% VTI = transversely isotropic medium with a vertical axis of symmetry

rho = 2500;  % density in kg/m^3
Vp = 3231;   % assumed isotropic medium P-wave velocity in m/s
Vs = 1844;   % assumed isotropic medium S-wave velocity in m/s
c11_0 = rho*Vp^2;
c44_0 = rho*Vs^2;
c12_0 = c11_0 - 2*c44_0;
C0 = [
   c11_0 c12_0 c12_0    0     0     0;
   c12_0 c11_0 c12_0    0     0     0;
   c12_0 c12_0 c11_0    0     0     0;
     0     0     0   c44_0    0     0;
     0     0     0      0   c44_0   0;
     0     0     0      0     0   c44_0;
    ];
 
% A-parameters for VTI media (Pšenčík and Farra, 2024)

    eps_x = 0.103;
    eps_z = 0;
    gamma_x = 0;
    gamma_z = 0.106;
    eta = -0.283;

% Elastic parameters of the VTI background
    c11_b = C0(1,1)*(1+2*eps_x);
    c12_b = C0(1,1)*(1+2*eps_x) - 2*C0(4,4)*(1+2*gamma_z);
    c13_b = C0(1,1)*(1+eps_x+eps_z+eta)- 2*C0(4,4)*(1+2*gamma_x);
    c33_b = C0(1,1)*(1+2*eps_z);
    c44_b = C0(4,4)*(1+2*gamma_x);
    c66_b = C0(4,4)*(1+2*gamma_z);

% Elastic parameters of the heterogeneous VFTI/orthorhombic media
    c11 = c11_b .* (1 - delN);
    c12 = c12_b .* (1 - delN);
    c13 = c13_b .* (1 - delN);
    c22 = c11_b .* (1 - delN .* (c12_b/c11_b)^2);
    c23 = c13_b .* (1 - delN .* (c12_b/c11_b));
    c33 = c33_B .* (1 - delN .* ((c13_b)^2 / (c33_b*c11_b)));
    c44 = c44_b .* ones(size(delN));   % the yz-plane is parallel to the fracture-plane
    c55 = c44_b .* (1 - delV);
    c66 = c66_b .* (1 - delH);
%% Case 1 ends here

%% Case 2: Heterogeneous VTI background (assuming we know the density and velocity models of the background):
% Comment all lines for Case 1 if you have background models else comment all lines for Case 2

load density_model.mat rho % density
load velocity_model.mat c11_b c33_b c44_b c66_b c13_b % VTI background stiffness coefficients
c12_b = c11_b - 2.*c66_b;

% Elastic parameters of the heterogeneous VFTI media
c11 = c11_b .* (1 - delN);
c12 = c12_b .* (1 - delN);
c13 = c13_b .* (1 - delN);
c22 = c11_b .* (1 - delN .* (c12_b ./ c11_b).^2);
c23 = c13_b .* (1 - delN .* (c12_b ./ c11_b));
c33 = c33_b .* (1 - delN .* ((c13_b.^2) ./ (c33_b .* c11_b)));
c44 = c44_b;
c55 = c44_b .* (1 - delV);
c66 = c66_b .* (1 - delH);
%% Case 2 ends here

save VFTI_model c11 c12 c13 c22 c23 c33 c44 c55 c66     % saves the model in a file for further use