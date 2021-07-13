%% Run a simple warm up problem
% Clean up and load libraries
clear all; clc; close all;

%% Path to mex function (c++ solver)
addpath(fullfile(pwd, "../../../build/lib"));

% Simple QP
nz = 2;
nx = 1;

Q = [2 0 0; 
     0 2 0;
     0 0 1];
g = [-2; -5; 0.1];

ubz = [1;6];

A = eye(nz+nx);
lbA = zeros(nz+nx, 1);
ubA = 100*ones(nz+nx, 1);

% Algorithm parameters (C++)
params.x0 = [1;];
params.z0 = [1; 3];
params.initialPenaltyParameter = 0.01;
params.penaltyUpdateFactor = 2;
params.qpSolver = 1;

% Run solver
miqp = struct('nz', nz, ...
              'nx', nx, ...
              'H', Q, ...
              'g', g, ...
              'ubz', ubz, ...
              'A', A, ...
              'lbA', lbA, ...
              'ubA', ubA);

 miqpsol(miqp, params);