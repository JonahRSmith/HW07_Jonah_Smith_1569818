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
close all;

%Part 1 - Done in Word
%Part 2 - Done below

TIMEO = 0; %Initial time
TIMEND = 10; %Final time, "T"
L = pi; %Length of 1-D 
D = 0.1; %Diffusivity
g0 = 0;
gL = 0;
k = 1;

Nodes = 100; %Control number of nodes across the domain 0 <= x <= L
dx = L/(Nodes-1);
x = dx:dx:L-dx;
TimeSteps = 1000; %Set the number of time steps your solution routine takes
dt = (TIMEND - TIMEO)/TimeSteps;

u = sin(k*x); %Define initial condition for u; u holds the current sol'n for u
uprev = u; %uprev holds the solution for u at the previous timestep


TIMEN = TIMEO; %Current time
CONST = D*dt/2/dx/dx; %Constant appearing many times through code
DNODES = length(x); %Number of diffusion nodes (i.e. not boundary nodes)
while (TIMEN < TIMEND)
    %Forward through the tri-diagonal!
    alpha = zeros(DNODES,1);
    g = zeros(DNODES,1);
    a = (1 + 2*CONST); %constant in all rows here
    b = -CONST*ones(DNODES,1); b(1)=0;
    c = -CONST*ones(DNODES,1); c(DNODES)=0;
    alpha(1) = a;
    g(1) = uprev(1) + CONST*(2*g0 - 2*uprev(1) + uprev(2));
    for j=2:DNODES-1
        alpha(j) = a - (b(j)*c(j-1))/alpha(j-1);
        g(j) = (uprev(j) + CONST*(uprev(j-1) - 2*uprev(j) + uprev(j+1))) - (b(j)*g(j-1)/alpha(j-1));
    end
    alpha(DNODES) = a - (b(j)*c(j-1))/alpha(j-1);
    g(DNODES) = (uprev(DNODES) + CONST*(uprev(DNODES-1) - 2*uprev(DNODES) + 2*gL)) - (b(j)*g(DNODES-1)/alpha(DNODES-1));
    u(DNODES) = g(DNODES)/alpha(DNODES);
    for kn=1:DNODES-1
        u(DNODES-kn) = (g(DNODES-kn) - c(DNODES-kn)*u(DNODES-kn+1))/alpha(DNODES-kn);
    end
    TIMEN = TIMEN + dt;
    uprev=u;
%     if (TIMEN+dt > TIMEND)
%         fprintf('Reducing dt from %g to %g for final iteration.\n',dt,TIMEND-TIMEN);
%         dt = TIMEND-TIMEN; %Make sure we end on TIMEND, and dont overshoot
%         CONST = D*dt/2/dx/dx; %Careful of makinig CONST too small
%     end
end

figure(1);
plot([0,x,L],[g0,u,gL]);
figure(2);
plot(x,exp(-D*k^2*TIMEN)*sin(k*x));

fprintf('END OF SCRIPT REACHED.\n');

