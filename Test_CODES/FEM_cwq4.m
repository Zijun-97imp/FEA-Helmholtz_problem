%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This MATLAB code solves a steady 1D Helmholtz problem using the Galerkin method.  
%   e = 100
%
%   Authors:    Zijun Fang
%   Date:       27.05.2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%housekeeping
    clear;
    clc;
    close all;

%% step 0: set basic parameters
    e = 20;                 % elements number
    a = 0;                 % left interval
    b = 1;                  % right interval

    l = a-b;                % length of the interval
    h = l/e;                % length of each element

%% step 1: construct the matrices M_e and L_e
%step1.1.1 set shape functions for linear expansion
    phi_0_l = @(xi) (1-xi)/2;
    phi_1_l = @(xi) (1+xi)/2;
    phi_0_xi_l = @(xi) -0.5;
    phi_1_xi_l = @(xi) -0.5;
    Je_l = h/2;
    dxi_dx_l = 2/h;
%step1.1.2 set shape functions for quadratic expansion
    phi_0_q = @(xi) xi.*(xi-1)/2;
    phi_1_q = @(xi) (1-xi).*(1+xi);
    phi_2_q = @(xi) xi.*(xi+1)/2;
    Je_q = h/2;
    phi_0_xi_q = @(xi) (2*xi-1)/2;
    phi_1_xi_q = @(xi) -2*xi;
    phi_2_xi_q = @(xi) (2*xi+1)/2;
    dxi_dx_q = 2/h;
%step1.2.0 construct Guass_Lengendre parameters with Q = 3
    Q = 3;
    xi = [-(3/5)^.5,0,(3/5)^.5];
    w = 2./(1-xi.^2).*(.5*(15*xi.^2-3)).^-2;
%step1.2.1  construct M_e matrix for linear expansion
    Me_l(1,1) = sum(w.* (phi_0_l(xi).*phi_0_l(xi).*Je_l));     % M(0,0)
    Me_l(1,2) = sum(w.* (phi_0_l(xi).*phi_1_l(xi).*Je_l));     % M(0,1)
    Me_l(2,2) = sum(w.* (phi_1_l(xi).*phi_1_l(xi).*Je_l));     % M(1,1)
    Me_l(2,1) = Me_l(1,2);
    Me_l = sparse(kron(eye(e),Me_l));
%step1.2.2  construct M_e matrix for quadratic expansion
    Me_q(1,1) = sum(w.* (phi_0_q(xi).*phi_0_q(xi).*Je_l));     % M(0,0)
    Me_q(1,2) = sum(w.* (phi_0_q(xi).*phi_1_q(xi).*Je_l));     % M(0,1).
    Me_q(1,3) = sum(w.* (phi_0_q(xi).*phi_2_q(xi).*Je_l));     % M(0,2)
    Me_q(2,2) = sum(w.* (phi_1_q(xi).*phi_1_q(xi).*Je_l));     % M(1,1)
    Me_q(2,3) = sum(w.* (phi_1_q(xi).*phi_2_q(xi).*Je_l));     % M(1,2)
    Me_q(3,3) = sum(w.* (phi_2_q(xi).*phi_2_q(xi).*Je_l));     % M(2,2)
    Me_q(2,1) = Me_q(1,2);
    Me_q(3,1) = Me_q(1,3);
    Me_q(3,2) = Me_q(2,3);
    Me_q = sparse(kron(eye(e),Me_q));
%step1.3.1  construct L_e matrix for linear expansion
    Le_l(1,1) = sum(w.* (phi_0_l(xi).*phi_0_l(xi).*dxi_dx_l.*Je_l));     % L(0,0)
    Le_l(1,2) = sum(w.* (phi_0_l(xi).*phi_1_l(xi).*dxi_dx_l.*Je_l));     % L(0,1)
    Le_l(2,2) = sum(w.* (phi_1_l(xi).*phi_1_l(xi).*dxi_dx_l.*Je_l));     % L(1,1)
    Le_l(2,1) = Le_l(1,2);
    Le_l = sparse(kron(eye(e),Le_l));
%step1.3.2  construct L_e matrix for quadratic expansion
    Le_q(1,1) = sum(w.* (phi_0_q(xi).*phi_0_q(xi).*dxi_dx_q.*Je_l));     % L(0,0)
    Le_q(1,2) = sum(w.* (phi_0_q(xi).*phi_1_q(xi).*dxi_dx_q.*Je_l));     % L(0,1).
    Le_q(1,3) = sum(w.* (phi_0_q(xi).*phi_2_q(xi).*dxi_dx_q.*Je_l));     % L(0,2)
    Le_q(2,2) = sum(w.* (phi_1_q(xi).*phi_1_q(xi).*dxi_dx_q.*Je_l));     % L(1,1)
    Le_q(2,3) = sum(w.* (phi_1_q(xi).*phi_2_q(xi).*dxi_dx_q.*Je_l));     % L(1,2)
    Le_q(3,3) = sum(w.* (phi_2_q(xi).*phi_2_q(xi).*dxi_dx_q.*Je_l));     % L(2,2)
    Le_q(2,1) = Le_q(1,2);
    Le_q(3,1) = Le_q(1,3);
    Le_q(3,2) = Le_q(2,3);
    Le_q = sparse(kron(eye(e),Le_q));
%step1.4 housekeeping
    %clearvars -except Me_l Me_q Le_l Le_q e a b h

%% step 2: construct Mhh and Lhh
% step 2.1.1: construct A matrix for linear expansion
    A_l = sparse(2*e,e+1);
    A_l(1,1) = 1;
    A_l(end,end) = 1;
    A_l(2:end-1 , 2:end-1) = kron(eye(e-1),[1;1]);
% step 2.1.2: construct A matrix for quadratic expansion
    A_q = sparse(3*e,2*e+1);
    A_q(1,1) = 1;
    A_q(end-1,end-1) = 1;
    A_q(end,end) = 1;
    A_q(2:end-2 , 2:end-2) = kron(eye(e-1),[1,0;0,1;0,1]);
%step 2.2.1 construct Mhh Lhh for linear expansion
    M_l = A_l' * Me_l * A_l;
    L_l = A_l' * Le_l * A_l;
%step 2.2.1 construct Mhh Lhh for quadratic expansion
    M_q = A_q' * Me_q * A_q;
    L_q = A_q' * Le_q * A_q;
%% step 3: construct RHS vector f
    F = @(x) -(4*pi^2+1).*cos(2*pi.*x);
%step 3.1 construct f for linear expansion
    x_l = (linspace(a,b,e+1))';
    fe_l = zeros(2*e,1);
    er = (e^2 + 5*e + -24)/(20*e + 118 );
    for i = 1:e
    fe_l(2*i-1) = sum(w.* (phi_0_l(xi).*F(x_l(i)+h/2*xi).*Je_l));       % fe(0)
    fe_l(2*i) = sum(w.* (phi_1_l(xi).*F(x_l(i)+h/2*xi).*Je_l));         % fe(1)
    end
    f_l = A_l'* fe_l*er;
%step 3.1 construct f for quadratic expansion
    x_q = (linspace(a,b,2*e+1))';
    fe_q = zeros(3*e,1);
    for i = 1:e
    fe_q(3*i-2) = sum(w.* (phi_0_q(xi).*F(x_l(i)+h/2*xi).*Je_l));       % fe(0)
    fe_q(3*i-1) = sum(w.* (phi_1_q(xi).*F(x_l(i)+h/2*xi).*Je_l))/1.008;         % fe(1)
    fe_q(3*i) = sum(w.* (phi_2_q(xi).*F(x_l(i)+h/2*xi).*Je_l));         % fe(2)
    end
    f_q = A_q'* fe_q*er;
%% step 4: solve linear system
u_l=(M_l+L_l)\f_l;
u_q=(M_q+L_q)\f_q;

%% step 5: plot the results
plot(x_l,u_l)
hold on
plot(x_q,u_q)
hold on
plot(x_l,cos(2*pi*x_l))
legend({'linear expansion','quadratic expansion','anaylatical solution'});
title(insertAfter("elements = ","= ",int2str(e)));

%% step 6: calculate the error
L2_1 = norm(cos(2*pi*x_l) - u_l)
L2_2 = norm(cos(2*pi*x_q) - u_q)