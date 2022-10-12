function [ke, me] = matelem_beam(E, I, L, S, rho, opt)

% Function OUTPUTS
% ke: matrice de raideur elementaire d'un element poutre
% me: matrice de masse d'un element poutre

% Function INPUTS
% E: module d'Young du materiau
% I: moment d'inertie de la section de la poutre
% L: Longueur de l'element
% S: aire de la section de la poutre 
% rho: masse volumique du matériau
% opt: 1 - matrice de masse consistante
%      2 - matrice de masse concentrée
%      3 - matrice de masse diagonale


C = E*I/L^3;

ke = C*[12 6*L -12 6*L; ...
       6*L 4*L^2 -6*L 2*L^2; ...
       -12 -6*L 12 -6*L; ...
       6*L 2*L^2 -6*L 4*L^2];
   
if opt == 1
    
    M = rho*S*L/420; 
    
    me = M*[156 22*L 54 -13*L; ...
           22*L 4*L^2 13*L -3*L^2; ...
           54 13*L 156 -22*L; ...
           -13*L -3*L^2 -22*L 4*L^2];
       
elseif opt == 2
      
    M = rho*S*L/2;
    
    me = M*diag([1 0 1 0]);
    
elseif opt == 3
    
    M = rho*S*L/78;
    
    me = M*diag([39 L^2 39 L^2]);
end