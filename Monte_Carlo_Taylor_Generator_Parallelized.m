% Mean Taylor Factor Test
%
% This script estimates mean Taylor factors and R-values by repeatedly
% sampling Euler orientations and minimizing the mean Taylor factor with
% respect to the lateral shrinking ratio q.
%
% Copyright (C) 2026 National Technology & Engineering Solutions of Sandia, LLC (NTESS). 
% Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
%
clc
clear
close all

% Read all orientations in the format Euler angles in radians (columns) for each grain (rows)
default_orientations = readmatrix('data.csv','FileType','text');

% Set up the parallel pool
numWorkers = 8; % Specify the number of workers (adjust based on your system)
if isempty(gcp('nocreate')) % Check if a parallel pool is already running
    parpool(numWorkers); % Create a parallel pool with the specified number of workers
end

%%
desiredIterations = 10000; % Number of desired iterations per number of grains.
grains_parameter = 1000*(1:5); % Grain counts explored.
for j = 1:length(grains_parameter)
    number_grains = grains_parameter(j);
    fprintf('Starting %.2f Grains\n', number_grains)
    csvFileName = strcat('Taylor_Model_',num2str(number_grains),'_Grains_Parallel.csv');
    
    parfor i = 1:desiredIterations
        selected_rows = datasample(default_orientations, number_grains, 'Replace', true); 
       
        fun = @(q) Calculate_Orientation(q,selected_rows);
        [q(i), M(i)] = fminsearch(fun,0.3); % Lateral Shrinking ratio to satisfy minimum Taylor Factor
        R(i) = q(i)/(1-q(i)); % Lankford Ratios (R-Value)
    
        if mod(i, round(desiredIterations * 0.05)) == 0
            fprintf('Progress: %.2f%%\n', (i / desiredIterations) * 100);
        end
    end
    
    iter = 1:desiredIterations;
    Results = [iter',M',R'];
    writematrix(Results,csvFileName)
end
delete(gcp('nocreate')); % Close the parallel pool

%%

function [M] = Calculate_Orientation(q,orientations)
    for i = 1:length(orientations)
        m(i) = FCC_Taylor_Factor(orientations(i,1),orientations(i,2),orientations(i,3),q,0.1);
    end
    M = mean(m);
end

