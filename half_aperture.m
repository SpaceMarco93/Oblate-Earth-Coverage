function [P1,P2,epsilon_1,epsilon_2] = half_aperture(a_tilde,b_tilde,alpha_SC,eta,r_line_local)

% half_aperture.m - Function to compute the area covered starting from the half-aperture angle.
%
% PROTOTYPE:
%   [P1,P2,epsilon_1,epsilon_2] = half_aperture(a_tilde,b_tilde,alpha_SC,eta,r_line_local)
%
% DESCRIPTION:
%   This function computes the instantaneous access area of a S/C for a 
%   given half-aperture angle.
%
% INPUT:
%   a_tilde             Semi-major axis of the intersected ellipse [km]
%   b_tilde             Semi-minor axis of the intersected ellipse [km]
%   alpha_SC            Angle of the S/C w.r.t. the axis e in local frame
%   eta                 Half-aperture angle [deg]
%   r_line_local        S/C position vector in the ellipse frame [km]
%
% OUTPUT:
%   P1                  Point of intersection on the right side [km]
%   P2                  Point of intersection on the left side [km]
%   epsilon_1           Elevation angle of the point on right side [deg]
%   epsilon_2           Elevation angle of the point on left side [deg]
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

% Definition of the S/C coordinates in the local frame
r_line_local = r_line_local(:);
e_sc = r_line_local(1);
u_sc = r_line_local(2);

% Aperture angles initialisation
eta = eta*pi/180;

% First case
if  alpha_SC >= 0 && alpha_SC <= pi/2  % First Quadrant
    
    % Angle of the secants w.r.to the semi-major axis direction
    alpha_P1 = (alpha_SC - eta);
    alpha_P2 = (alpha_SC + eta);
    
    % Slopes of the two secants
    m_P1 = tan(alpha_P1);
    m_P2 = tan(alpha_P2);
    
    % Coefficients to normalise the equation and reduce the error
    coeff2_1 = (2*m_P1*a_tilde^2*(u_sc - m_P1*e_sc))/(a_tilde^2*m_P1^2 + b_tilde^2);
    coeff3_1 = (-a_tilde^2*b_tilde^2 + a_tilde^2*(u_sc-m_P1*e_sc)^2)/(a_tilde^2*m_P1^2 + b_tilde^2);
    
    coeff2_2 = (2*m_P2*a_tilde^2*(u_sc - m_P2*e_sc))/(a_tilde^2*m_P2^2 + b_tilde^2);
    coeff3_2 = (-a_tilde^2*b_tilde^2 + a_tilde^2*(u_sc-m_P2*e_sc)^2)/(a_tilde^2*m_P2^2 + b_tilde^2);
    
    % Discriminant of the second-degree equation
    Delta_P1 = coeff2_1^2 - 4*coeff3_1;
    Delta_P2 = coeff2_2^2 - 4*coeff3_2;
    
    % Selection of the two right solutions among the possible four
    if alpha_P2 <= pi/2
        
        e_P1 = (-coeff2_1 + sqrt(Delta_P1))/2;
        e_P2 = (-coeff2_2 + sqrt(Delta_P2))/2;
        u_P1 = m_P1*e_P1 - m_P1*e_sc + u_sc;
        u_P2 = m_P2*e_P2 - m_P2*e_sc + u_sc;
        
        % Points in vector form
        P1 = [e_P1,u_P1,0];
        P2 = [e_P2,u_P2,0];
        
    end
    if alpha_P2 > pi/2
        
        e_P1 = (-coeff2_1 + sqrt(Delta_P1))/2;
        e_P2 = (-coeff2_2 - sqrt(Delta_P2))/2;
        u_P1 = m_P1*e_P1 - m_P1*e_sc + u_sc;
        u_P2 = m_P2*e_P2 - m_P2*e_sc + u_sc;
        
        % Points in vector form
        P1 = [e_P1,u_P1,0];
        P2 = [e_P2,u_P2,0];
        
    end
end

% Second case
if  alpha_SC > pi/2 && alpha_SC <= pi   % Second Quadrant
    
    % Angle of the secants w.r.to the semi-major axis direction
    alpha_P1 = (alpha_SC - eta);
    alpha_P2 = (alpha_SC + eta);
    
    % Slopes of the two secants
    m_P1 = tan(alpha_P1);
    m_P2 = tan(alpha_P2);
    
    % Coefficients to normalise the equation and reduce the error
    coeff2_1 = (2*m_P1*a_tilde^2*(u_sc - m_P1*e_sc))/(a_tilde^2*m_P1^2 + b_tilde^2);
    coeff3_1 = (-a_tilde^2*b_tilde^2 + a_tilde^2*(u_sc-m_P1*e_sc)^2)/(a_tilde^2*m_P1^2 + b_tilde^2);
    
    coeff2_2 = (2*m_P2*a_tilde^2*(u_sc - m_P2*e_sc))/(a_tilde^2*m_P2^2 + b_tilde^2);
    coeff3_2 = (-a_tilde^2*b_tilde^2 + a_tilde^2*(u_sc-m_P2*e_sc)^2)/(a_tilde^2*m_P2^2 + b_tilde^2);
    
    % Discriminant of the second-degree equation
    Delta_P1 = coeff2_1^2 - 4*coeff3_1;
    Delta_P2 = coeff2_2^2 - 4*coeff3_2;
    
    % Selection of the two right solutions among the possible four
    if alpha_P1 <= pi/2
        
        e_P1 = (-coeff2_1 + sqrt(Delta_P1))/2;
        e_P2 = (-coeff2_2 - sqrt(Delta_P2))/2;
        u_P1 = m_P1*e_P1 - m_P1*e_sc + u_sc;
        u_P2 = m_P2*e_P2 - m_P2*e_sc + u_sc;
        
        % Points in vector form
        P1 = [e_P1,u_P1,0];
        P2 = [e_P2,u_P2,0];
        
    end
    if alpha_P1 > pi/2
        
        e_P1 = (-coeff2_1 - sqrt(Delta_P1))/2;
        e_P2 = (-coeff2_2 - sqrt(Delta_P2))/2;
        u_P1 = m_P1*e_P1 - m_P1*e_sc + u_sc;
        u_P2 = m_P2*e_P2 - m_P2*e_sc + u_sc;
        
        % Points in vector form
        P1 = [e_P1,u_P1,0];
        P2 = [e_P2,u_P2,0];
        
    end
end

% Third case
if  alpha_SC >= -pi && alpha_SC <= -pi/2   % Third Quadrant
    
    % Angle of the secants w.r.to the semi-major axis direction
    alpha_P1 = (alpha_SC - eta);
    alpha_P2 = (alpha_SC + eta);
    
    % Slopes of the two secants
    m_P1 = tan(alpha_P1);
    m_P2 = tan(alpha_P2);
    
    % Coefficients to normalise the equation and reduce the error
    coeff2_1 = (2*m_P1*a_tilde^2*(u_sc - m_P1*e_sc))/(a_tilde^2*m_P1^2 + b_tilde^2);
    coeff3_1 = (-a_tilde^2*b_tilde^2 + a_tilde^2*(u_sc-m_P1*e_sc)^2)/(a_tilde^2*m_P1^2 + b_tilde^2);
    
    coeff2_2 = (2*m_P2*a_tilde^2*(u_sc - m_P2*e_sc))/(a_tilde^2*m_P2^2 + b_tilde^2);
    coeff3_2 = (-a_tilde^2*b_tilde^2 + a_tilde^2*(u_sc-m_P2*e_sc)^2)/(a_tilde^2*m_P2^2 + b_tilde^2);
    
    % Discriminant of the second-degree equation
    Delta_P1 = coeff2_1^2 - 4*coeff3_1;
    Delta_P2 = coeff2_2^2 - 4*coeff3_2;
    
    % Selection of the two right solutions among the possible four
    if alpha_P2 <= -pi/2
        
        e_P1 = (-coeff2_1 - sqrt(Delta_P1))/2;
        e_P2 = (-coeff2_2 - sqrt(Delta_P2))/2;
        u_P1 = m_P1*e_P1 - m_P1*e_sc + u_sc;
        u_P2 = m_P2*e_P2 - m_P2*e_sc + u_sc;
        
        % Points in vector form
        P1 = [e_P1,u_P1,0];
        P2 = [e_P2,u_P2,0];
        
    end
    if alpha_P2 > -pi/2
        
        e_P1 = (-coeff2_1 - sqrt(Delta_P1))/2;
        e_P2 = (-coeff2_2 + sqrt(Delta_P2))/2;
        u_P1 = m_P1*e_P1 - m_P1*e_sc + u_sc;
        u_P2 = m_P2*e_P2 - m_P2*e_sc + u_sc;
        
        % Points in vector form
        P1 = [e_P1,u_P1,0];
        P2 = [e_P2,u_P2,0];
        
    end
end

% Fourth case
if  alpha_SC > -pi/2 && alpha_SC < 0     % Fourth quadrant
    
    % Angle of the secants w.r.to the semi-major axis direction
    alpha_P1 = (alpha_SC - eta);
    alpha_P2 = (alpha_SC + eta);
    
    % Slopes of the two secants
    m_P1 = tan(alpha_P1);
    m_P2 = tan(alpha_P2);
    
    % Coefficients to normalise the equation and reduce the error
    coeff2_1 = (2*m_P1*a_tilde^2*(u_sc - m_P1*e_sc))/(a_tilde^2*m_P1^2 + b_tilde^2);
    coeff3_1 = (-a_tilde^2*b_tilde^2 + a_tilde^2*(u_sc-m_P1*e_sc)^2)/(a_tilde^2*m_P1^2 + b_tilde^2);
    
    coeff2_2 = (2*m_P2*a_tilde^2*(u_sc - m_P2*e_sc))/(a_tilde^2*m_P2^2 + b_tilde^2);
    coeff3_2 = (-a_tilde^2*b_tilde^2 + a_tilde^2*(u_sc-m_P2*e_sc)^2)/(a_tilde^2*m_P2^2 + b_tilde^2);
    
    % Discriminant of the second-degree equation
    Delta_P1 = coeff2_1^2 - 4*coeff3_1;
    Delta_P2 = coeff2_2^2 - 4*coeff3_2;
    
    % Selection of the two right solutions among the possible four
    if alpha_P1 <= -pi/2
        
        e_P1 = (-coeff2_1 - sqrt(Delta_P1))/2;
        e_P2 = (-coeff2_2 + sqrt(Delta_P2))/2;
        u_P1 = m_P1*e_P1 - m_P1*e_sc + u_sc;
        u_P2 = m_P2*e_P2 - m_P2*e_sc + u_sc;
        
        % Points in vector form
        P1 = [e_P1,u_P1,0];
        P2 = [e_P2,u_P2,0];
        
    end
    if alpha_P1 > -pi/2
        
        e_P1 = (-coeff2_1 + sqrt(Delta_P1))/2;
        e_P2 = (-coeff2_2 + sqrt(Delta_P2))/2;
        u_P1 = m_P1*e_P1 - m_P1*e_sc + u_sc;
        u_P2 = m_P2*e_P2 - m_P2*e_sc + u_sc;
        
        % Points in vector form
        P1 = [e_P1,u_P1,0];
        P2 = [e_P2,u_P2,0];
        
    end
end

%% Determination of the tangents in P1 and P2

% Tangent slope
if u_P1 >= 0           % Selection of the semi-ellipse
    m_t_P1 = (-b_tilde*e_P1)/(a_tilde^2*sqrt(1 - (e_P1/a_tilde)^2));
else
    m_t_P1 = (b_tilde*e_P1)/(a_tilde^2*sqrt(1 - (e_P1/a_tilde)^2));
end

if u_P2 >= 0
    m_t_P2 = (-b_tilde*e_P2)/(a_tilde^2*sqrt(1 - (e_P2/a_tilde)^2));
else
    m_t_P2 = (b_tilde*e_P2)/(a_tilde^2*sqrt(1 - (e_P2/a_tilde)^2));
end

% Elevation angle computation: angle between two lines
epsilon_1 = atan((-m_t_P1 + m_P1)/(1 + m_t_P1*m_P1));
epsilon_2 = atan((m_t_P2 - m_P2)/(1 + m_t_P2*m_P2));

% Conversion into degrees
epsilon_1 = epsilon_1*180/pi;
epsilon_2 = epsilon_2*180/pi;

end