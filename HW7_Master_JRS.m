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
figure(1);
ct=0;
while (TIMEN < TIMEND)
    ct = ct+1;
    %Forward through the tri-diagonal!
    alpha = zeros(DNODES);
    g = zeros(DNODES);
    a = (1 + 2*CONST); %constant in all rows here
    b = -CONST; %Constant, except where undefined in row(1)
    c = -CONST; %Constant, except where undefined in row(DNODES)
    alpha(1) = 1 + 2*CONST;
    g(1) = uprev(1) + CONST*(2*g0 - 2*uprev(1) + uprev(2));
    for j=2:DNODES-1
        alpha(j) = a - (b*c)/alpha(j-1);
        g(j) = (uprev(j) + CONST*(uprev(j-1) - 2*uprev(j) + uprev(j+1))) - (b*g(j-1)/alpha(j-1));
    end
    alpha(DNODES) = a - (b*c)/alpha(j-1);
    g(DNODES) = uprev(DNODES) + CONST*(uprev(DNODES-1) - 2*uprev(DNODES) + 2*gL);
    u(DNODES) = g(DNODES)/alpha(DNODES);
    for k=1:DNODES-1
        u(DNODES-k) = (g(DNODES-k) - c*u(DNODES-k))/alpha(DNODES-k);
    end
    TIMEN = TIMEN + dt;
    uprev=u;
    if (TIMEN+dt > TIMEND)
        fprintf('Reducing dt from %g to %g for final iteration.\n',dt,TIMEND-TIMEN);
        dt = TIMEND-TIMEN; %Make sure we end on TIMEND, and dont overshoot
    end
    if (ct == 10)
        ct = 0;
        hold on;
        plot([0,x,L],[g0,u,gL]);
    end
end
fprintf('END OF SCRIPT REACHED.\n');






