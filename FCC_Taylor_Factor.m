function [M] = FCC_Taylor_Factor(phi_1,Phi, phi_2, q, eta)
% FCC_TAYLOR_FACTOR Compute the Taylor factor for FCC slip systems.
%
% Copyright (C) MTEX authors/contributors.
% Copyright (C) 2026 National Technology & Engineering Solutions of Sandia, LLC (NTESS). 
% Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
%
% Adapted from MTEX function calcTaylor.m (MTEX Toolbox, https://mtex-toolbox.github.io/)
% for standalone FCC Taylor-factor calculation on 2026-09-08.
%
% Licensed under GNU GPL v2 or later; see LICENSE in repository root.

    %% Strain Transformation
    E_s = eta *[1, 0, 0; 0, -q, 0; 0, 0, -(1-q)]; % Strain Tensor in Specimen Coordinates
    
    a = zeros(3,3);
    a(1,1) = cos(phi_1)*cos(phi_2) - sin(phi_1)*sin(phi_2)*cos(Phi);
    a(1,2) = cos(phi_2)*sin(phi_1) + sin(phi_2)*cos(phi_1)*cos(Phi);
    a(1,3) = sin(phi_2)*sin(Phi);
    a(2,1) = -cos(phi_1)*sin(phi_2) - sin(phi_1)*cos(phi_2)*cos(Phi);
    a(2,2) = -sin(phi_2)*sin(phi_1) + cos(phi_1)*cos(phi_2)*cos(Phi);
    a(2,3) = cos(phi_2)*sin(Phi);
    a(3,1) = sin(phi_1)*sin(Phi);
    a(3,2) = -cos(phi_1)*sin(Phi);
    a(3,3) = cos(Phi);
    
    eps = zeros(3, 3); % Strain Tensor in Crystal Coordinates
    for i = 1:3
        for j = 1:3
            eps(i,j) = a(i,1) * a(j,1) * E_s(1,1) + ...
                        a(i,1) * a(j,2) * E_s(1,2) + ...
                        a(i,1) * a(j,3) * E_s(1,3) + ...
                        a(i,2) * a(j,1) * E_s(2,1) + ...
                        a(i,2) * a(j,2) * E_s(2,2) + ...
                        a(i,2) * a(j,3) * E_s(2,3) + ...
                        a(i,3) * a(j,1) * E_s(3,1) + ...
                        a(i,3) * a(j,2) * E_s(3,2) + ...
                        a(i,3) * a(j,3) * E_s(3,3);
        end
    end
    
    %% Coefficients for FCC slip systems
    e = zeros(3,3,24);
    e(1,1,:) = [ 0 -1  1  0  1 -1  1  1  0 -1 -1  0  0 -1 -1  0  1  1 -1  1  0  1 -1  0];
    e(1,2,:) = [ 1 -1  0 -1  1  0  0  1  1  0 -1 -1 -1  1  0  1 -1  0  1  0  1 -1  0 -1]*0.5;
    e(1,3,:) = [-1  0  1  1  0 -1 -1  0  1  1  0 -1  1  0  1 -1  0 -1  0  1  1  0 -1 -1]*0.5;
    e(2,1,:) = [ 1 -1  0 -1  1  0  0  1  1  0 -1 -1 -1  1  0  1 -1  0  1  0  1 -1  0 -1]*0.5;
    e(2,2,:) = [ 1  0 -1 -1  0  1 -1  0  1  1  0 -1  1  0  1 -1  0 -1  0 -1 -1  0  1  1];
    e(2,3,:) = [ 0  1 -1  0 -1  1  1  1  0 -1 -1  0  0  1  1  0 -1 -1 -1  1  0  1 -1  0]*0.5;
    e(3,1,:) = [-1  0  1  1  0 -1 -1  0  1  1  0 -1  1  0  1 -1  0 -1  0  1  1  0 -1 -1]*0.5;
    e(3,2,:) = [ 0  1 -1  0 -1  1  1  1  0 -1 -1  0  0  1  1  0 -1 -1 -1  1  0  1 -1  0]*0.5;
    e(3,3,:) = [-1  1  0  1 -1  0  0 -1 -1  0  1  1 -1  1  0  1 -1  0  1  0  1 -1  0 -1];        
    e = e/sqrt(6);
    
    %% Determine the slip on each system.
    % initialize the coefficients
    b = zeros(1,length(e));
    
    % critical resolved shear stress - CRSS
    % by now assumed to be identical - might also be stored in sS
    CRSS = ones(length(e),1);
    
    % decompose eps into sum of disclocation tensors, that is we look for
    % coefficients b such that sSepsSym * b = eps
    
    % since the strain tensor is symmetric we require only 5 entries out of it
    A = reshape(e,9,[]);
    A = A([1,2,3,5,6],:);
    
    % the strain coefficients to match
    y = reshape(eps,9,[]);
    y = y([1,2,3,5,6],:);
    
    options = optimoptions('linprog','Algorithm','interior-point','Display','none');
    
    % for all strain tensors do
    for i = 1:size(y,2)
      
      % determine coefficients b with A * b = y and such that sum |CRSS_j *
      % b_j| is minimal. This is equivalent to the requirement b>=0 and CRSS*b
      % -> min which is the linear programming problem solved below
          b(i,:) = linprog(CRSS,[],[],A,y(:,i),zeros(size(A,2),1),[],options);      
    end
    
    %% Calculating the Taylor Factor
    % the Taylor factor is simply the sum of the coefficents
    M = sum(b,2) ./ norm(eps);
end