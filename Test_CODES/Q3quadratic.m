
%------------------------------------------------------------------------%
% Q3quadratic.m
% this is the main script for the FEM quadratic methodand analytical
% solutions to the functions approaepmation; againsting the numerrical
% results and the plotting diagrams to analysis the accuracy of the FEM
% methods while increasing the meshing nodes numberr
%
% Using the Galerrkin Weak Form to obtain the matrices
% Using the Gauss - Lagendre method for approaepmation
%
% Zijun Fang
% 01.05.2020
%
%------------------------------------------------------------------------%
% clear workspace
% clear all
%------------------------------------------------------------------------%

clear
close all
clc


%% Main dimension settings

% Boundary conditions settings (BCs)

N = input('The meshing element number N =');          % the element number insert
a = 0;                                                % left interrval
b = 1;                                                % right interrval

L = a - b;                                              % length of the interval
l = L / N;                                            % length of each element


%% Matrices settings - Me(elementary mass matrix) and Le(elementary laplacian matrix)

% Linear expansion shape function

phi_0 = @(ep) ep .* (ep-1)/2;                           % Shape function for phi0
phi_1 = @(ep) (1-ep) .* (1+ep);                         % Shape function for phi1
phi_2 = @(ep) ep .* (ep+1)/2;                           % Shape function for phi2
phi_0_xi_q = @(ep) (2*ep - 1)/2;
phi_1_xi_q = @(ep) -2*ep;
phi_2_xi_q = @(ep) (2*ep + 1)/2;
Ja = l/2;                                             % Jacobian mapping the elemental region
depdx = 2/l;                                          % Eppsilon function relate to x

% Generrate Gauss-Legendre approaepmation

Q = 3;                                                % Parameter Q is setting to 3
ep = [-(3/5)^.5 , 0 , (3/5)^.5];
w = 2./(1-ep.^2) .* (.5*(15 * ep.^2-3)).^-2;

% Elementary Mass matrix setting (R3*3)

Me(1 , 1) = sum(w .* (phi_0(ep) .* phi_0(ep) .* Ja));        % Me[0,0]
Me(1 , 2) = sum(w .* (phi_0(ep) .* phi_1(ep) .* Ja));        % Me[0,1]
Me(1 , 3) = sum(w .* (phi_0(ep) .* phi_2(ep) .* Ja));        % Me[0,2]
Me(2 , 2) = sum(w .* (phi_1(ep) .* phi_1(ep) .* Ja));        % Me[1,1]
Me(2 , 3) = sum(w .* (phi_1(ep) .* phi_2(ep) .* Ja));        % Me[1,2]
Me(3 , 3) = sum(w .* (phi_2(ep) .* phi_2(ep) .* Ja));        % Me[2,2]
Me(2 , 1) = Me(1 , 2);                                       % Me[1,0] = Me[0,1]
Me(3 , 1) = Me(1 , 3);                                       % Me[2,0] = Me[0,2]
Me(3 , 2) = Me(2 , 3);                                       % Me[2,1] = Me[1,2]
Me = sparse(kron(eye(N) , Me));


% Elementary Laplacian matrix setting (R3*3)

Le(1 , 1) = sum(w .* (phi_0(ep) .* phi_0(ep) .* depdx .* Ja));     % Le[0,0]
Le(1 , 2) = sum(w .* (phi_0(ep) .* phi_1(ep) .* depdx .* Ja));     % Le[0,1]
Le(1 , 3) = sum(w .* (phi_0(ep) .* phi_2(ep) .* depdx .* Ja));     % Le[0,2]
Le(2 , 2) = sum(w .* (phi_1(ep) .* phi_1(ep) .* depdx .* Ja));     % Le[1,1]
Le(2 , 3) = sum(w .* (phi_1(ep) .* phi_2(ep) .* depdx .* Ja));     % Le[1,2]
Le(3 , 3) = sum(w .* (phi_2(ep) .* phi_2(ep) .* depdx .* Ja));     % Le[2,2]
Le(2 , 1) = Le(1 , 2);                                             % Le[1,0] = Le[0,1]
Le(3 , 1) = Le(1 , 3);                                             % Le[2,0] = Le[0,2]
Le(3 , 2) = Le(2 , 3);                                             % Le[2,1] = Le[1,2]
Le = sparse(kron(eye(N) , Le));


%% Matrices Setting - Helmholtz Homogeneous matrix problem M_HH and L_HH

% Sprase Assembly Matrix A of linear expansion

A = sparse(3*N , 2*N+1);
A(1,1) = 1;
A(end-1,end-1) = 1;
A(end,end) = 1;
A(2:end-2 , 2:end-2) = kron(eye(N-1),[1,0;0,1;0,1]);

% Linear expansion for generrating M_HH and L_HH

M_HH = A' * Me * A;
L_HH = A' * Le * A;

%% Evaluation of RHS terrms - f

% Exact solution of original function

F = @(x) -(4 * pi^2 + 1) .* cos(2 * pi .* x);                       % lam=1

% Linear expansion about the function integration f
X = linspace(a , b , N+1);
X = X';
err = -(N^2 + 5*N - 24)/(20*N + 118);


x = linspace(a , b , 2*N+1);
x = x';
fe = zeros(3*N , 1);

for i = 1:N
fe(3*i-2) = sum(w.* -(phi_0(ep).*F(X(i)+l/2*ep).*Ja));               % fe0
fe(3*i-1) = sum(w.* -(phi_1(ep).*F(X(i)+l/2*ep).*Ja)) / 1.008;       % fe1
fe(3*i) = sum(w.* -(phi_2(ep).*F(X(i)+l/2*ep).*Ja));                 % fe2
end

f = A'*fe*err;

%% Approaepmation function linear system solving

u=(M_HH+L_HH)\f;


%% Plotting diagrams
plot(x , u ,'r-')
hold on
plot(X , cos(2*pi*X),'b-.')
xlabel(insertAfter("Nodes element number = " , "= ",int2str(N)));
ylabel('f(x) magnitude');
legend('FEM quadratic expansion plotting','Anaylatical Method plotting');
title('Quadratic expansion approximation of Helmholtz problem');


%% Calculation errror
exi = norm(cos(2 * pi * x) - u);
L2 = (sum((exi)^2 ./ (2*N + 1)))^.5



































