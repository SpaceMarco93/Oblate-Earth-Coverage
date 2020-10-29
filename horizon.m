function [eta_hor_1,eta_hor_2,lambda_hor,alpha_SC] = horizon(a_tilde,b_tilde,r_line_local)

% horizon.m - Function to compute the horizon coverage properties.
%
% PROTOTYPE:
%   [eta_hor_1,eta_hor_2,lambda_hor,alpha_SC] = horizon(a_tilde,b_tilde,r_line_local)
%
% DESCRIPTION:
%   This function computes the horizon coverage properties that are the
%   quantities obtained considering an ideal scenario without the
%   foreshortening effect. The properties are computed considering the
%   tangents to the ellipse starting from the S/C position.
%
% INPUT:
%   a_tilde             Semi-major axis of the intersected ellipse [km]
%   b_tilde             Semi-minor axis of the intersected ellipse [km]
%   r_line_local        S/C position vector in the ellipse frame [km]
%
% OUTPUT:
%   eta_hor_1           Horizon-boresight angle of the right side [rad]
%   eta_hor_2           Horizon-boresight angle of the left side [rad]
%   lambda_hor          Horizon ground-rage angle [rad]
%   alpha_SC            Angle of the S/C w.r.t. the axis e in local frame
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

% Definition of the S/C coordinates in the local frame
r_line_local = r_line_local(:);
e_sc = r_line_local(1);
u_sc = r_line_local(2);

% Coefficients derived from the normalisation of the second-degree
% equation to reduce te numerical error
coeff_1 = (2*e_sc*u_sc)/(a_tilde^2 - e_sc^2);
coeff_2 = (b_tilde^2 - u_sc^2)/(a_tilde^2 - e_sc^2);

% Discrimant of the second degree-equation
Delta_hor = coeff_1^2 - 4*coeff_2;

% Slopes of the tangents: condition to assign the proper P1 and P2
if  e_sc*u_sc >= 0
    m_T1 = (-coeff_1 + sqrt(Delta_hor))/2;
    m_T2 = (-coeff_1 - sqrt(Delta_hor))/2;
else
    m_T1 = (-coeff_1 - sqrt(Delta_hor))/2;
    m_T2 = (-coeff_1 + sqrt(Delta_hor))/2;
end

% Vertical intercept of the tangents
q_T1 = -m_T1*e_sc + u_sc;
q_T2 = -m_T2*e_sc + u_sc;

% Horizon Points
e_T1 = (m_T1*(-q_T1)*a_tilde^2)/(b_tilde^2 + a_tilde^2*m_T1^2);
e_T2 = (m_T2*(-q_T2)*a_tilde^2)/(b_tilde^2 + a_tilde^2*m_T2^2);
u_T1 = m_T1*e_T1 + q_T1;
u_T2 = m_T2*e_T2 + q_T2;
r_T1 = [e_T1;u_T1;0];
r_T2 = [e_T2;u_T2;0];

% Assign the position vector P1 and P2 to the forward and backward points
cross_prod = cross(r_T2,r_T1);
if cross_prod(3) < 0
    r_T1 = [e_T2;u_T2;0]; 
    r_T2 = [e_T1;u_T1;0];
end

% Boresight-angle and ground range-angle
eta_hor_1 = acos(((r_line_local-r_T1)'*r_line_local)/(norm(r_line_local-r_T1)*norm(r_line_local)));
eta_hor_2 = acos(((r_line_local-r_T2)'*r_line_local)/(norm(r_line_local-r_T2)*norm(r_line_local)));
lambda_hor = acos(r_T1'*r_T2/(norm(r_T1)*norm(r_T2)));

% Angle between the S/C line-of-sight and the the semi-major axis direction e
alpha_SC = atan2(u_sc,e_sc);

end