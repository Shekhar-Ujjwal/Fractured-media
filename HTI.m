%  Code to generate stiffness coefficients for the heterogeneous HTI (transversely isotropic media with horizontal axis of symmetry)
%  model based on the linear slip theory (Schoenberg, 1983)
% Author: Ujjwal Shekhar

% Input fracture weaknesses at grid points 
load fracture_weakness_model.mat delN delT

% delN = normal fracture weakness, delT = tangential fracture weakness
% Fracture weaknesses usually varies from 0 to 0.5 in rock mass, should not be  >=1 or < 0

[Nx, Ny, Nz] = size(delN);  % Extracting the size of the model (can be both 3D and 2D)
 
%% Case 1: Homogeneous isotropic background (in which vertical fractures are embedded):
rho = 2500;  % density in kg/m^3
Vp = 3231;   % P-wave velocity in m/s
Vs = 1844;   % S-wave velocity in m/s
c11_b = rho*Vp^2;
c44_b = rho*Vs^2;
c12_b = c11_b - 2*c44_b;
C0 = [
   c11_b c12_b c12_b    0     0     0;
   c12_b c11_b c12_b    0     0     0;
   c12_b c12_b c11_b    0     0     0;
     0     0     0   c44_b    0     0;
     0     0     0      0   c44_b   0;
     0     0     0      0     0   c44_b;
    ];
 
% Elastic parameters of the heterogeneous HTI media
c11 = C0(1,1) .* (1 - delN);
c12 = C0(1,2) .* (1 - delN);
c13 = C0(1,3) .* (1 - delN);
c22 = C0(1,1) .* (1 - delN .* (C0(1,2)/C0(1,1))^2);
c23 = C0(1,3) .* (1 - delN .* (C0(1,2)/C0(1,1)));
c33 = C0(3,3) .* (1 - delN .* ((C0(1,3))^2 / (C0(3,3)*C0(1,1))));
c44 = C0(4,4) .* ones(size(delN));   % as the yz-plane is set parallel to planes of vertical fractures
c55 = C0(4,4) .* (1 - delT);
c66 = C0(6,6) .* (1 - delT);
%% Case 1 ends here

%% Case 2: Heterogeneous isotropic background (assuming we know the density and velocity models of the background):
% Comment all lines for Case 1 if you have background models else comment all lines for Case 2

load density_model.mat rho % density
load Vp_model.mat Vp       % P-wave velocity
load Vs_model.mat Vs       % S-wave velocity
c11_b = rho.*Vp^2;
c44_b = rho.*Vs^2;
c12_b = c11_b - 2.*c44_b;
c13_b = c12_b;

% Elastic parameters of the heterogeneous HTI media
c11 = c11_b .* (1 - delN);
c12 = c12_b .* (1 - delN);
c13 = c13_b .* (1 - delN);
c22 = c11_b .* (1 - delN .* (c12_b ./ c11_b).^2);
c23 = c13_b .* (1 - delN .* (c12_b ./ c11_b));
c33 = c33_b .* (1 - delN .* ((c13_b.^2) ./ (c33_b .* c11_b)));
c44 = c44_b;
c55 = c44_b .* (1 - delT);
c66 = c44_b .* (1 - delT);
%% Case 2 ends here

save HTI_model c11 c12 c13 c22 c23 c33 c44 c55 c66     % saves the model in a file for further use