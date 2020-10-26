function [r_SC_proj,eta_hor,lambda_hor,P1_in,P2_in] = coverage_function(r_SC,angle,n,flag,N,tol)

% coverage_function.m - Instantaneous access area in the Geocentric frame for a generic perturbed pointing.
%
% PROTOTYPE:
%   [nadir,eta_hor,lambda_hor,P1_in,P2_in] = coverage_function(r_SC,angle,n,flag,N,tol)
%
% DESCRIPTION:
%   This function computes the instantaneous access area of a S/C at a given
%   height considering the Earth as an oblate ellipsoid of rotation and
%   using as input the half-aperture angle or the elevation angle for a
%   generic perturbed pointing.
%
% INPUT:
%   r_SC                S/C position vector in Geocentric Inertial frame [km]
%   angle               Starting angle:
%                           - half-aperture angle [deg]
%                           - elevation angle [deg]
%   n                   Line-of-sight direction in Geocentric frame
%   flag                Flag to decide the type of input:
%                           - 1 for half-aperture angle (default)
%                           - 0 for elevation angle
%   N                   Discretisation of the conical field-of-view
%   tol                 Tolerance for the elevation angle cycle:
%                           - default: 1e-4
%
% OUTPUT:
%   r_SC_proj           Projectorion of the S/C position on the Earth
%                       surface [km]
%   eta_hor             Horizon boresight angle [deg]
%   lambda_hor          Hprzon ground-range angle [deg]
%   P1_in               1st set of instantaneous access area points [km]
%   P2_in               2nd set of instantaneous access area points [km]
%
% FUNCTIONS CALLED:
%   nadir.m, ellipse_prop.m, horizon.m, half_aperture.m, elevation.m
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

%% Initialisation

% Assign the default parameters
if nargin < 6
    tol = 1e-4;
    if nargin < 5
        N = 30;
        if nargin < 4
            flag = 1;
        end
    end
end

% Control to check the type of input angle
if flag == 1
    eta = angle;
else
    epsilon = angle;
end

% Geometric data for the Earth oblate ellipsoid (WGS-84)
a = 6378.1363;                        % Oblate ellipsoid semi-major axis [km]
b = 6356.7516005;                     % Oblate ellipsoid semi-minor axis [km]
E = sqrt(1 - b^2/a^2);                % Oblate ellipsoid eccentricity

% Fix the S/C position vector and the direction as column vectors
r_SC = r_SC(:);
n = n(:);

%% Line-of-sight projection onto the Earth surface (nadir)
r_SC_proj = nadir(r_SC,a,b,n);

%% Definition of a new reference frame
% The new reference frame is shifted w.r.t. the Geocentric frame:
% - origin in the intersection of the line of sight with the Equatorial
%   Plane
% - axes parallel to the Geocentric Inertial frame

% Opposite direction to the line of sight
o = -n;

% Determination of the geographic coordinates
lon_nadir = atan2(o(2),o(1));
lat_nadir = asin(o(3));

% Definition of the angles of the rotation matrix 321
phi = lon_nadir;
theta = -lat_nadir;
psi = linspace(0,pi,N);

% Initialisation of the varibles
P1_in = zeros(N,3);
P2_in = zeros(N,3);

% Rotation of the generic plane around the line of sight
for i = 1:N
    
    % Rotation matrix from Geocentric to horizon frame
    A_321 = [cos(theta)*cos(phi),                                   cos(theta)*sin(phi),                                   -sin(theta);
        -cos(psi(i))*sin(phi)+sin(psi(i))*sin(theta)*cos(phi),  cos(psi(i))*cos(phi)+sin(psi(i))*sin(theta)*sin(phi),   sin(psi(i))*cos(theta);
        sin(psi(i))*sin(phi)+sin(theta)*cos(phi)*cos(psi(i)), -sin(psi(i))*cos(phi) + sin(theta)*sin(phi)*cos(psi(i)), cos(theta)*cos(psi(i))];
    
    % Normal vector to the generic plane
    n_prime = A_321'*[0;0;1];
    
    % Distance of the plane w.r.to the origin
    d = n_prime'*r_SC;
    
    % Geometrical properties of the intersected ellipse
    [a_tilde,b_tilde,s,~,~,A_rot,r_line_local] = ellipse_prop(r_SC,a,b,E,n_prime,d);
    
    %% Horizon determination: tangent points
    
    % Horizon ground-range and half-aperture angles
    [eta_hor_1,eta_hor_2,lambda_hor,alpha_SC] = horizon(a_tilde,b_tilde,r_line_local);
    
    % Total horizon boresight angle
    eta_hor = eta_hor_1 + eta_hor_2;
    
    %% Computation of the effective horizon
    
    % Compute the points in the local frame according to the method
    if flag == 1
        [P1,P2,~,~] = half_aperture(a_tilde,b_tilde,alpha_SC,eta,r_line_local);
    else
        [P1,P2,~,~] = elevation(a_tilde,b_tilde,alpha_SC,eta_hor_1,eta_hor_2,epsilon,tol,r_line_local);
    end
    
    % Inverse transformation from local to Geocentric inertial frame
    P1_3D = A_rot'*P1';
    P2_3D = A_rot'*P2';
    
    % Relative reference systems formulation to align the origin
    P1_in(i,:) = P1_3D + s;
    P2_in(i,:) = P2_3D + s;
    
end

% Conversion of angles in degrees
eta_hor = eta_hor*180/pi;
lambda_hor = lambda_hor*180/pi;

end
