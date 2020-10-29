function test()

% test.m  Script for the application of the coverage function.
%
% DESCRIPTION:
%   Script for the application of the coverage function determination.
%   The example is used considering as reference case a Galileo satellite.
%   The user can change the type of S/C setting different orbital paramenters
%   and replacing the default "Orbit" in the script with the equivalent orbit
%   obtained from the conversion of the orbital parameters into Cartesian
%   coordinates.
%   The coverage can be done using as input the half-aperture angle of the
%   conical field of view or the minimum elevation angle associated to the
%   locations on the Earth surfaces.
%
% FUNCTIONS CALLED:
%   coverage_function.m
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

close all

% WGS-84 parameters definition
a_Earth = 6378.137;                 % Earth semi-major axid [km]
f = 1/298.257223563;                % Earth flattening
mu_Earth = 398600.4418;             % Earth gravitational parameter [km^3/s^2]
w_Earth = 7.292115e-5;              % Earth rotation rate  WGS-84[rad/s]

% Derived quantities
R_eq = a_Earth;                     % Equatorial radius [km]
R_pol = a_Earth*(1-f);              % Polar radius [km]        
% R_pol = R_eq;                      % Polar radius [km]

% Define the coverage properties
eta = 10;                           % Half-aperture angle [deg]
epsilon = 10;                       % Minimum elevation angle [deg]
N = 50;                             % Discretisation of the conical signal

% Keplerian elements of a Galileo satellite
a = 29600;                          % Semi-major axis [km]
e = 0;                              % Eccentricity
i = 56*pi/180;                      % Inclination [rad]
Om = 240*pi/180;                    % Longitude of ascending node [rad]
om = 0;                             % Pericenter anomaly [rad]
ni = 0;                             % Initial true anomaly [rad]

% Orbital period
T = 2*pi*sqrt(a^3/mu_Earth);

% Define the vector of true anomalies
theta = linspace(ni,2*pi + ni,360);

% Compute the vector of eccentric anomalies
E = 2*atan(sqrt((1-e)/(1+e)).*tan(theta./2));

% Compute the time vector related to the true anomaly vector
t = sqrt(a^3/mu_Earth).*(E-e.*sin(E));
t(t < 0) = t(t < 0) + T;

% Construction of the Earth surface
n_panels = 180;
[x_Earth,y_Earth,z_Earth] = ellipsoid(0,0,0,R_eq, R_eq, R_pol, n_panels);

% Inertial coordinates of the Earth oblate ellipsoid
Earth_in = [x_Earth(:), y_Earth(:), z_Earth(:)];

% Load Earth image for texture map
Earth_image = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';
cdata = imread(Earth_image);

% Initialisation of variables
Earth_rot = zeros(length(x_Earth(:)),3);

% Load the Galileo orbit
load('Example.mat','Orbit');

% Plot
for j = 1:360
    
    % Define an indefinite plot
    h = plot3(nan,nan,nan,'or'); hold on;
    
    % Plot the S/C orbit
    plot3(Orbit(1,:),Orbit(2,:),Orbit(3,:),'Color',[0 0.4470 0.7410]);
    set(h,'XData',Orbit(1,j),'YData',Orbit(2,j),'ZData',Orbit(3,j))
    
    % Rotation Matrix to simulate the Earth rotation
    A_rot = [cos(w_Earth*t(j)), -sin(w_Earth*t(j)), 0;
        sin(w_Earth*t(j)), cos(w_Earth*t(j)), 0;
        0,                  0,                 1];
    
    for k = 1:size(Earth_in,1)
        % Rotation of the coordinates
        Earth_rot(k,:) = A_rot*Earth_in(k,:)';
    end
    
    % Define the S/C position vector at time j
    r_SC = Orbit(:,j);
    
    % Define the Geocentric direction at time j
    n_GEO = -r_SC/norm(r_SC);
    
    % Compute the coverage with one of the two methods
%    [~,~,~,P1_in,P2_in] = coverage_function(R_eq,R_pol,r_SC,eta,n_GEO,1,N);
    [~,~,~,P1_in,P2_in] = coverage_function(R_eq,R_pol,r_SC,epsilon,n_GEO,0,N);
    
    % Save the result of the Earth point in view of the satellite
    area = [P1_in;
        P2_in];
    
    % Plot the area covered by the S/C
    plot3(area(:,1),area(:,2),area(:,3),'-g');
    
    % Define the Earth rotation
    x_Earth_rot = reshape(Earth_rot(:,1),sqrt(size(Earth_rot(:,1),1)),sqrt(size(Earth_rot(:,1),1)));
    y_Earth_rot = reshape(Earth_rot(:,2),sqrt(size(Earth_rot(:,2),1)),sqrt(size(Earth_rot(:,2),1)));
    z_Earth_rot = reshape(Earth_rot(:,3),sqrt(size(Earth_rot(:,3),1)),sqrt(size(Earth_rot(:,3),1)));
    
    % Earth texture
    globe = surf(x_Earth_rot,y_Earth_rot,-z_Earth_rot,'FaceColor','none','EdgeColor',0.5*[1 1 1]);
    set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
    hold on;
    grid on;
    xlabel('X [km]');
    ylabel('Y [km]');
    zlabel('Z [km]');
    axis equal
    drawnow
    hold off
end



