%% HOMEWORK #07, Jonah R. Smith, 1569818
% Diffusion Equation and Crank-Nicolson Methods
%One dimensional problem
% du/dt = D * d^2u/dx^2
% 0 <= t <= T
% 0 <= x <= L
% Dirichlet boundary conditions:
% u(0,t) = g0; u(L,t) = gL
% Initial Condition:
% u(x,0) = f(x)

clc;
clear all;

%Part 1 - Done in Word
%Part 2 - Done below

TIMEO = 0; %Initial time
TIMEND = 10; %Final time, "T"
L = pi; %Length of 1-D 
D = 0.1; %Diffusivity
g0 = 0;
gL = 0;
k = 1;
f = @(x) sin(k*x);
