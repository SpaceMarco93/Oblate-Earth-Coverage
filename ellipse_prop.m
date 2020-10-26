function [a_tilde,b_tilde,s,e,u,A_rot,r_line_local] = ellipse_prop(r_SC,a,b,E,n_prime,d)

% ellipse_prop.m - Function to compute the geometric quantities of the
% local ellipse.
%
% PROTOTYPE:
%   [a_tilde,b_tilde,s,e,u,A_rot,r_line_local] = ellipse_prop(r_SC,a,b,E,n_prime,d)
%
% DESCRIPTION:
%   This function computes the geometric properties of the ellipse obtained
%   by the intersection of a generic plane with the Earth oblate ellipsoid
%   of rotation.
%
% INPUT:
%   r_SC                S/C position vector in Geocentri frame [km]
%   a                   Semi-major axis of the oblate ellipsoid [km]
%   b                   Semi-minor axis of the oblate ellipsoid [km]
%   E                   Eccentricity of the oblate ellipsoid
%   n_prime             Normal to the ellipse plane in the Geocentric frame
%   d                   Distance between the Earth centre and the ellipse
%                       plane [km]
%
% OUTPUT:
%   a_tilde             Semi-major axis of the intersected ellipse [km]
%   b_tilde             Semi-minor axis of the intersected ellipse [km]
%   s                   Position vector of the ellipse centre
%   e                   Direction of the semi-major axis
%   u                   Direction of the semi-minor axis
%   A_rot               Rotation matrix from the Geocentric Inertial frame
%                       to the local ellipse frame
%   r_line_local        S/C position vector in the ellipse frame
%
% FUNCTIONS CALLED:
%
% AUTHOR:
%   Marco Nugnes, 24/10/2020, https://www.compass.polimi.it,
%   E-mail: marco.nugnes@polimi.it
%      
% REFERENCE AND LICENSE: 
%   Copyright 2020 Marco Nugnes
%   https://www.compass.polimi.it
%
%   This set of codes is distributed under the 3-clause BSD license (see 
%   below) with the additional clause to cite the reference paper where the
%   theoretical work is explained and the website of the COMPASS project, 
%   which funded the research:
%   - Nugnes M., Colombo, C., and Tipaldi, M., "Coverage Area Determination
%	for Conical Fields of View Considering an Oblate Earth", Journal of
%	Guidance, Control, and Dynamics, Vol. 42, No. 10, pp. 2233-2245, 2019.
%	DOI: https://doi.org/10.2514/1.G004156.
%   - https://compass.polimi.it.
%
% ACKNWOLEDGEMENT
%   The research leading to these results has received funding from the 
%   European Research Council (ERC) under the European Unions Horizon 2020 
%   research and innovation program as part of project COMPASS 
%   (Grant agreement No. 679086)
%
% -----------------------------------------------------------------------

% Set the components of the normal to the ellipse plane
n1_prime = n_prime(1);
n2_prime = n_prime(2);
n3_prime = n_prime(3);

% Ellipse's centre in the Geocentric inertial frame
s = [(n1_prime*d)/(1-E^2*n3_prime^2), (n2_prime*d)/(1-E^2*n3_prime^2), ((b^2/a^2)*n3_prime*d)/(1-E^2*n3_prime^2)]';

% Ellipse's semi-major axis
a_tilde = a*sqrt(1 - d^2/(a^2*(1-E^2*n3_prime^2)));

% Ellipse's semi-minor axis
b_tilde = b*(sqrt(1 - (d^2/a^2) - E^2*n3_prime^2))/(1-E^2*n3_prime^2);

if n1_prime == 0 && n2_prime == 0 && n3_prime == 1
    
    % Ellipse's semi-major axis direction in the Geocentric frame
    e = [1,0,0];
    
    % Ellipse's semi-major axis direction in the Geocentric frame
    u = [0,1,0];
    
elseif n1_prime == 0 && n2_prime == 0 && n3_prime == -1
    
    % Ellipse's semi-major axis direction in the Geocentric frame
    e = [1,0,0];
    
    % Ellipse's semi-major axis direction in the Geocentric frame
    u = [0,-1,0];
    
else
    
    % Ellipse's semi-major axis direction in the Geocentric frame
    e = (1/sqrt(n1_prime^2 + n2_prime^2))*[n2_prime,-n1_prime,0];
    
    % Ellipse's semi-minor axis direction in the Geocentric frame
    u = -(1/sqrt(n1_prime^2+n2_prime^2))*[-n1_prime*n3_prime,-n2_prime*n3_prime,n1_prime^2+n2_prime^2];
    
end

% Position vector of the S/C w.r.to the ellipse's centre
r_line = r_SC - s;

% Rotation matrix from Geocentric to the Local frame of the ellipse
A_rot = [e;
    u;
    n_prime'];

% Position vector of the S/C w.r.to the ellipse centre
r_line_local = A_rot*r_line;

end