function [x, y] = thin(run, n, kthin)
% PURPOSE: function called by raftery.m

% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jpl@jpl.econ.utoledo.edu

% NOTE: this is a translation of FORTRAN code
%       (which is why it looks so bad)

 ind = 1:kthin:n; 
 y = run(1, ind);
 t1 = size(y, 2);
 x = t1;