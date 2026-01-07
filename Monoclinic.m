% Effective monoclinic model for two non-orthogonal vertical fracture sets in VTI background
% Author: Ujjwal Shekhar

delN1 = 0.10; delV1 = 0.15; delH1 = 0.20; % fracture weaknesses set 1
delN2 = 0.06; delV2 = 0.10; delH2 = 0.15; % fracture weaknesses set 2

phi1 = pi/4; phi2 = pi/2;   % azimuth angles of fracture sets

% Thomsen parameters based stiffness matrix elements of the background VTI medium %

rho = 2300;  % density in kg/m^3
Vp0 = 2600; Vs0 = 1200; % (P and S-waves vertical velocities in m/s)
C33b = rho*Vp0^2/1e9; C44b = rho*Vs0^2/1e9; % stiffnesses in GPa
epsilon = 0.1; gamma = 0.15; delta = 0.05; 

C11b = (2*epsilon +1)*C33b;
C66b = (2*gamma +1)*C44b;
C13b = sqrt(2*C33b*(C33b-C44b)*delta + (C33b-C44b)^2) -C44b;
C12b = C11b - 2*C66b;

% VTI background stiffness matrix
Cb = [C11b C12b C13b 0    0     0
      C12b C11b C13b 0    0     0
      C13b C13b C33b 0    0     0
      0    0    0    C44b 0     0
      0    0    0    0    C44b  0
      0    0    0    0    0     C66b];

% Fracture compliance matrices

Kn1 = delN1/(C11b*(1-delN1));  Kv1 = delV1/(C44b*(1-delV1));  Kh1 = delH1/(C66b*(1-delH1));  
Kn2 = delN2/(C11b*(1-delN2));  Kv2 = delV2/(C44b*(1-delV2));  Kh2 = delH2/(C66b*(1-delH2)); 

D1 = [Kn1  0    0    0    0     0
      0    0    0    0    0     0
      0    0    0    0    0     0
      0    0    0    0    0     0
      0    0    0    0    Kv1   0
      0    0    0    0    0     Kh1];

D2 = [Kn2  0    0    0    0     0
      0    0    0    0    0     0
      0    0    0    0    0     0
      0    0    0    0    0     0
      0    0    0    0    Kv2   0
      0    0    0    0    0     Kh2];

% Bond rotation matrices

R1 = [(cos(phi1))^2       (sin(phi1))^2     0    0          0          sin(2*phi1)
      (sin(phi1))^2       (cos(phi1))^2     0    0          0         -sin(2*phi1)
      0                   0                 1    0          0          0
      0                   0                 0    cos(phi1) -sin(phi1)  0
      0                   0                 0    sin(phi1)  cos(phi1)  0
      -0.5*sin(2*phi1)    0.5*sin(2*phi1)   0    0          0          cos(2*phi1)];

R2 = [(cos(phi2))^2       (sin(phi2))^2     0    0          0          sin(2*phi2)
      (sin(phi2))^2       (cos(phi2))^2     0    0          0         -sin(2*phi2)
      0                   0                 1    0          0          0
      0                   0                 0    cos(phi2) -sin(phi2)  0
      0                   0                 0    sin(phi2)  cos(phi2)  0
      -0.5*sin(2*phi2)    0.5*sin(2*phi2)   0    0          0          cos(2*phi2)];

% Stiffness matrix of monoclinic media with the horizontal symmetry plane
Cmono = Cb*inv(eye(6) + (R1*D1*R1' + R2*D2*R2')*Cb);  % the unit of elements of the matrix is GPa

% Density-normalized stiffness coefficients (useful for ray-based modelling methods)
% (Note that thin fractures have no significant effect on the density of rock )
A11 = Cmono(1,1)/rho; A12 = Cmono(1,2)/rho; A13 = Cmono(1,3)/rho; A16 = Cmono(1,6)/rho;
A22 = Cmono(2,2)/rho; A23 = Cmono(2,3)/rho; A26 = Cmono(2,6)/rho;
A33 = Cmono(3,3)/rho; A36 = Cmono(3,6)/rho;
A44 = Cmono(4,4)/rho; A45 = Cmono(4,5)/rho;
A55 = Cmono(5,5)/rho;
A66 = Cmono(6,6)/rho;
