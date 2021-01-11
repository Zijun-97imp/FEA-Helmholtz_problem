%------------------------------------------------------------------------%
% Q3linear.m
% this is the main script for the FEM lienar methodand analytical 
% solutions to the functions approaximation; againsting the numerical 
% results and the plotting diagrams to analysis the accuracy of the FEM 
% methods while increasing the meshing nodes number
% 
% Using the Galerkin Weak Form to obtain the matrices 
% Using the Gauss - Lagendre method for approaximation
%
% Zijun Fang
% 01.05.2020
%
%------------------------------------------------------------------------%
% clear workspace
% clear all
%------------------------------------------------------------------------%


clc
clear

%% Main dimension settings

% Boundary conditions settings (BCs)

a = 0;           % Left side interval
b = 1;           % Right side interval
N = input('Meshing element number N =');           % Meshing element number

% Axial dimension settings

L = abs(a - b);                % Simulation length
l = L / N;                % Mesh element length


%% Matrices settings - Me(elementary mass matrix) and Le(elementary laplacian matrix)

% Linear expansion shape function

phi_0 = @(ep) (1 - ep) / 2;                          % Shape function for phi0
phi_1 = @(ep) (1 + ep) / 2;                          % Shape function for phi1
Ja = l / 2;                                          % Jacobian mapping the elemental region
depdx = 2 / l;                                       % Eppsilon function relate to x

% Generate Gauss-Legendre approaximation

Q = 3;                                            % Parameter Q is setting to 3
ep = [-(3/5)^.5 , 0 , (3/5)^.5];
w = 2 ./ (1 - ep.^2) .* (.5 * (15 * ep.^2-3)).^-2;

% Elementary Mass matrix setting (R2*2)

Me(1 , 1) = sum(w .* (phi_0(ep) .* phi_0(ep) .* Ja));     % Me[0,0]
Me(1 , 2) = sum(w .* (phi_0(ep) .* phi_1(ep) .* Ja));     % Me[0,1]
Me(2 , 2) = sum(w .* (phi_1(ep) .* phi_1(ep) .* Ja));     % Me[1,1]
Me(2 , 1) = Me(1 , 2);                                    % Me[1,0] = Me[0,1]
Me = sparse(kron(eye(N) , Me));

% Elementary Laplacian matrix setting (R2*2)

Le(1 , 1) = sum(w .* (phi_0(ep) .* phi_0(ep) .* depdx .* Ja));     % Le[0,0]
Le(1 , 2) = sum(w .* (phi_0(ep) .* phi_1(ep) .* depdx .* Ja));     % Le[0,1]
Le(2 , 2) = sum(w .* (phi_1(ep) .* phi_1(ep) .* depdx .* Ja));     % Le[1,1]
Le(2 , 1) = Le(1 , 2);                                               % Le[1,0] = Le[0,1]
Le = sparse(kron(eye(N) , Le));


%% Matrices Setting - Helmholtz Homogeneous matrix problem M_HH and L_HH

% Sprase Assembly Matrix A of linear expansion

A = sparse(2*N , N+1);                                 % Sprase matrix A (non-orthogonal) size
A(1 , 1) = 1;
A(end , end) = 1;
A(2:end-1 , 2:end-1) = kron(eye(N-1) , [1;1]);         % Linear expansion matrix setting

% Linear expansion for generating M_HH and L_HH

M_HH = A' * Me * A;              % M = A'*Me*A
L_HH = A' * Le * A;              % L = A'*Le*A



%% Evaluation of RHS terms - f

% Exact solution of original function

F = @(x) -(4 * pi^2 + 1) * cos(2 * pi * x);                           % lamda = 1

% Linear expansion about the function integration f

x = linspace(a , b , N+1);
x = x';
fe = zeros(2*N , 1);
err = -(N^2 + 5*N + -24)/(20*N + 118);

for i = 1 : N
    fe(2*i-1) = sum(w .* (phi_0(ep) .* F(x(i) + l/2*ep) .* Ja));       % fe0
    fe(2*i) = sum(w .* (phi_1(ep) .* F(x(i) + l/2*ep) .* Ja));         % fe1
end
f = A' * fe * err;


%% Approaximation function linear system solving

u = (M_HH + L_HH) \ f;

%% Plotting diagrams
plot(x , u,'r-')
hold on
plot(x , cos(2*pi*x),'b.-')
xlabel(insertAfter("Nodes element number = " , "= ",int2str(N)));
ylabel('f(x) magnitude');
legend('FEM linear expansion plotting','Anaylatical Method plotting');
title('Linear expansion approximation of Helmholtz problem');


%% Calculation error
exi = norm(cos(2 * pi * x) - u);
L2 = (sum((exi)^2 ./ (N + 1)))^.5









