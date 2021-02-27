# FEA-Helmholtz_problem
The individual project for FEA method of Helmholtz problem.

+ %-------------------------------------------%
  %README.m
  %
  %this is the code running folder, including
  %the figure output file 'Figure'
  %and the code source file 'Code'. 
  %
  %
  %Zijun Fang
  %01.05.2020
  %
  %-------------------------------------------%


+ MATLAB base running file

  % FEM linear and quadratic expansion coursework

  % figure output-loglog plot

  % figure output-node against plot



The concentrated euqation of Helmholtz Problem:
                         (d/dx)*( (phy)*(du/dx) ) - (lam)u = f(x);
                         
Neumann and Dirichlet boundary conditions:
                             u(a) = alpha; du/dx(b) = beta;

Alternative mixed boundary condition:
                           (phy) [(du/dx)(a)] + u(a) = gamma;

* Number of mesh: N_e = 5, 10, 20, 50, 100;
