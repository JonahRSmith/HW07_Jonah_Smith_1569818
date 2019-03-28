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

%% 2.a

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
dt = (TIMEND - TIMEO)/TimeSteps; nomdt=dt;

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
    if (TIMEN+dt > TIMEND)
        if (TIMEN < TIMEND)
            fprintf('Reducing dt from %g to %g for final iteration(s).\n',dt,TIMEND-TIMEN);
            dt = TIMEND-TIMEN; %Make sure we end on TIMEND, and dont overshoot
            CONST = D*dt/2/dx/dx; %Careful of makinig CONST too small
        end
    end
end

figure(1);
hold on;
title(sprintf('Part (a), k=%i, N=%i, dt=%g, TIMEN=%0.3f',k,DNODES+2,nomdt,TIMEN));
plot([0,x,L],[g0,u,gL], 'b-');
exact_u = exp(-D*k^2*TIMEN)*sin(k*x);
plot([0,x,L],[g0,exact_u,gL], 'r--');
legend('Crank-Nicolson Solution', 'Real Solution');

%Calculating average error:
err = 0;
for ni = 1:DNODES
    err = err + abs((u(ni) - exact_u(ni))/exact_u(ni));
end
err = err/DNODES;
text(0.125,max(u),sprintf('Average Error: %g',err));

fprintf('END OF SCRIPT REACHED.\n');

%% 2.b

clear all;

TIMEO = 0; %Initial time
TIMEND = 10; %Final time, "T"
L = pi; %Length of 1-D 
D = 0.1; %Diffusivity
%g0 = sin(omega*t);
%gL = sin(omega*t)*cos(k*L);
%F(x,t) = (omega*cos(omega*t) + D*k^2*sin(omega*t))*cos(k*x);
omega = 1;
k = 1;

Nodes = 100; %Control number of nodes across the domain 0 <= x <= L
dx = L/(Nodes-1);
x = dx:dx:L-dx;
TimeSteps = 1000; %Set the number of time steps your solution routine takes
dt = (TIMEND - TIMEO)/TimeSteps; nomdt=dt;

u = 0*x; %Define initial condition for u; u holds the current sol'n for u
uprev = u; %uprev holds the solution for u at the previous timestep


TIMEN = TIMEO; %Current time
CONST = D*dt/2/dx/dx; %Constant appearing many times through code
DNODES = length(x); %Number of diffusion nodes (i.e. not boundary nodes)
while (TIMEN < TIMEND)
    %Forward through the tri-diagonal!
    g0 = sin(omega*TIMEN);
    gL = sin(omega*TIMEN)*cos(k*L);
    alpha = zeros(DNODES,1);
    g = zeros(DNODES,1);
    a = (1 + 2*CONST); %constant in all rows here
    b = -CONST*ones(DNODES,1); b(1)=0;
    c = -CONST*ones(DNODES,1); c(DNODES)=0;
    alpha(1) = a;
    g(1) = uprev(1) + CONST*(2*g0 - 2*uprev(1) + uprev(2));
    for j=2:DNODES-1
        alpha(j) = a - (b(j)*c(j-1))/alpha(j-1);
        F = (omega*cos(omega*TIMEN) + D*k^2*sin(omega*TIMEN))*cos(k*x(j));
        g(j) = dt*F + (uprev(j) + CONST*(uprev(j-1) - 2*uprev(j) + uprev(j+1))) - (b(j)*g(j-1)/alpha(j-1));
    end
    alpha(DNODES) = a - (b(j)*c(j-1))/alpha(j-1);
    g(DNODES) = (uprev(DNODES) + CONST*(uprev(DNODES-1) - 2*uprev(DNODES) + 2*gL)) - (b(j)*g(DNODES-1)/alpha(DNODES-1));
    u(DNODES) = g(DNODES)/alpha(DNODES);
    for kn=1:DNODES-1
        u(DNODES-kn) = (g(DNODES-kn) - c(DNODES-kn)*u(DNODES-kn+1))/alpha(DNODES-kn);
    end
    TIMEN = TIMEN + dt;
    uprev=u;
    if (TIMEN+dt > TIMEND)
        if (TIMEN < TIMEND)
            fprintf('Reducing dt from %g to %g for final iteration(s).\n',dt,TIMEND-TIMEN);
            dt = TIMEND-TIMEN; %Make sure we end on TIMEND, and dont overshoot
            CONST = D*dt/2/dx/dx; %Careful of makinig CONST too small
        end
    end
end

figure(2);
hold on;
title(sprintf('Part (a), k=%i, N=%i, dt=%g, TIMEN=%0.3f, omega=%g',k,DNODES+2,nomdt,TIMEN,omega));
plot([0,x,L],[g0,u,gL], 'b-');
exact_u = sin(omega*TIMEN)*cos(k*x);
plot([0,x,L],[g0,exact_u,gL], 'r--');
legend('Crank-Nicolson Solution', 'Real Solution');

%Calculating average error:
err = 0;
for ni = 1:DNODES
    err = err + abs((u(ni) - exact_u(ni))/exact_u(ni));
end
err = err/DNODES;
text(0.125,max(u),sprintf('Average Error: %g',err));

fprintf('END OF SCRIPT REACHED.\n');

