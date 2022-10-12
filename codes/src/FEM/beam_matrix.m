function [Ks, Ms, bmesh] = beam_matrix(xmin, L, S, I, E, rho, Nelem, opt)

% [Ks, Ms, bmesh] = beam_matrix(xmin, L, E, I,S, rho) construit les matrices de
% masse et de raideur d'une poutre
%
% Entrees:
%    - xmin: origine du maillage FEM [m]
%    - L: Longueur de la poutre [m]
%    - S: Aire de la section droite [m^2]
%    - I: Moment quadratique de la section [m^4]
%    - E: Module d'Young du materiau [Pa]
%    - rho: Masse volumique du materiau [kg/m^3]
%    - Nelem: Nombre d'element du maillage
%    - opt: Definie le type de CL appliquees aux extremites 
%           * opt = 1: simplement appuyee - simplement appuye
%           * opt = 2: libre - libre
%           * opt = 3: encastre - encastre
%
% Sorties:
%    - Ks: Matrice de raideur de la poutre 
%    - Ms: Marice de masse de la poutre
%    - bmesh: Maillage de la poutre

bmesh = beam_mesh(xmin, L, Nelem);                            % Maillage de la poutre
[ke, me] = matelem_beam(E, I, bmesh.elem_size, S, rho, 1); % Construction des matrices elementaires
K = matrix_assembly(ke, bmesh);                            % Construction de la matrice de raideur
M = matrix_assembly(me, bmesh);                            % Construction de la matrice de masse

Nnodes = size(bmesh.Nodes, 1);

if opt == 1
    bmesh.nCL = [2:2*Nnodes-2, 2*Nnodes];  % ddls non contraints
    
elseif opt == 2
    bmesh.nCL = 1:2*Nnodes;
    
elseif opt == 3
    bmesh.nCL = 3:2*Nnodes-2;
end

Ks = K(bmesh.nCL, bmesh.nCL);
Ms = M(bmesh.nCL, bmesh.nCL);