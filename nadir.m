function r_SC_proj = nadir(r_SC,a,b,n)

% nadir.m - Function to compute the projection of the S/C on the Earth surface.
%
% PROTOTYPE:
%   r_SC_proj = nadir(r_SC,a,b,n)
%
% DESCRIPTION:
%   This function is used to compute the projection of the S/C on the Earth
%   surface along the navigation signal direction.
%
% INPUT:
%   r_SC                S/C position vector in Geocentri frame [km]
%   a                   Semi-major axis of the oblate ellipsoid [km]
%   b                   Semi-minor axis of the oblate ellipsoid [km]
%   n                   Pointing direction of the navigation signal in the
%                       Geocentric inertial frame
%
% OUTPUT:
%   r_SC_proj           Projection of the S/C center of mass on the Earth
%                       surface
%
% FUNCTIONS CALLED:
%
% AUTHOR:
%   Marco Nugnes, 24/10/2020, https://www.compass.polimi.it,
%   E-mail: marco.nugnes@polimi.it
%      
% REFERENCE AND LICENSE: 
%   Copyright 2020 Marco Nugnes
%   This code is made available under the Creative Commons 
%   Attribution-NonCommercial-ShareAlike 4.0 International(CC BY-NC-SA 4.0)
%   This license is accessible at:
%   https://creativecommons.org/licenses/by-nc-sa/4.0/
%   The code is free to use for research purposes, but whenever used I 
%   kindly ask to cite the following article where the theoretical 
%   framework of the code is explained:
%   Nugnes M., Colombo, C., and Tipaldi, M., "Coverage Area Determination 
%   for Conical Fields of View Considering an Oblate Earth", Journal of 
%   Guidance, Control, and Dynamics, Vol. 42, No. 10, pp. 2233-2245, 2019.
%   DOI: https://doi.org/10.2514/1.G004156.
%   For more info about this research visit the website: 
%   https://compass.polimi.it. 
%   For commercial use, please contact the author. 
%
% ACKNWOLEDGEMENT
%   The research leading to these results has received funding from the 
%   European Research Council (ERC) under the European Unions Horizon 2020 
%   research and innovation program as part of project COMPASS 
%   (Grant agreement No. 679086)
%
% -----------------------------------------------------------------------

% Coordinates of the S/C in the Geocentric Inertial Frame
XX = r_SC(1);
YY = r_SC(2);
ZZ = r_SC(3);

% Direction cosines of the line of sight
nx = n(1);
ny = n(2);
nz = n(3);

% Discriminant of the second-degree equation
Delta_tan = 2*b^4*nx*ny*XX*YY+2*a^2*b^2*nx*nz*XX*ZZ+2*a^2*b^2*ny*nz*YY*ZZ...
    -b^4*nx^2*YY^2-a^2*b^2*nx^2*ZZ^2 + a^2*b^4*nx^2 - b^4*ny^2*XX^2 ...
    -a^2*b^2*ny^2*ZZ^2 + a^2*b^4*ny^2 - a^2*b^2*nz^2*XX^2 ...
    -a^2*b^2*YY^2*nz^2 + a^4*b^2*nz^2;

% Solution of the parametric equation
t_1 = (-(b^2*nx*XX+b^2*ny*YY+a^2*nz*ZZ) + sqrt(Delta_tan))/(b^2*nx^2 + b^2*ny^2 + a^2*nz^2);
t_2 = (-(b^2*nx*XX+b^2*ny*YY+a^2*nz*ZZ) - sqrt(Delta_tan))/(b^2*nx^2 + b^2*ny^2 + a^2*nz^2);

% Selection of the minimum value corresponding to the real projection
t = min(t_1,t_2);

% Computation of the Intersection point coordinates
Int_x = r_SC(1) + nx*t;
Int_y = r_SC(2) + ny*t;
Int_z = r_SC(3) + nz*t;
r_SC_proj = [Int_x,Int_y,Int_z]';

end