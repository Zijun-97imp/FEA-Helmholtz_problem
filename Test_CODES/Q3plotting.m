%------------------------------------------------------------------------%
% Q3plotting.m
% 
% this script is running the log-log diagram to plot the error for both the
% linear and quadratic finite element expansions, versus the size of the
% element (1/Nel);
%
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
close all



%% General dimensional settings
% Meshing element size (1/Nel)

Nel = [5 10 20 50 100];                        % Input element number
s = Nel.^-1;                                   % Size of the element

% FEM error norm number L2

L2_l = [0.4928 0.2469 0.1210 0.0462 0.0222];   % Norm number of the error in linear expansion
L2_q = [0.6543 0.3133 0.1272 0.0447 0.0229];   % Norm number of the error in quadratic expansion

% log-log function

S = log(s);
ER_l = log(L2_l);
ER_q = log(L2_q);


%% Plotting log-log diagram

% log-log diragram plotting

figure(1)
plot(S , ER_l , 'r+-' , 'markersize' , 10)
hold on
plot(S , ER_q , 'bx-.', 'markersize' , 10)
xlabel('Size of meshing element size - log(1/Nel)', 'FontSize' , 14);
ylabel('FEM analysis expansion error - log(L2)', 'FontSize' , 14);
legend('FEM linear expansion error' , 'FEM quadratic expansion error', 'FontSize' , 12);
title('FEM Helmholtz prblem analysis linear and quadratic expansion error log-log plotting', 'FontSize' , 15);
grid on

figure(2)
plot(s , L2_l , 'r+-' , 'markersize' , 10)
hold on
plot(s , L2_q , 'bx-.', 'markersize' , 10)
xlabel('Size of meshing element size - 1/Nel' , 'FontSize' , 14);
ylabel('FEM analysis expansion error - L2', 'FontSize' , 14);
legend('FEM linear expansion error' , 'FEM quadratic expansion error', 'FontSize' , 12);
title('FEM Helmholtz prblem analysis linear and quadratic expansion error plotting', 'FontSize' , 15);
grid on





















