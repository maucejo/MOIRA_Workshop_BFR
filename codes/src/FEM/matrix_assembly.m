function K = matrix_assembly(ke, smesh)

% Function OUTPUTS
% K: matrice assemblée

% Function INPUTS
% ke: matrice elementaire d'un element
% smesh: maillage de la structure

Nnodes = size(smesh.Nodes, 1);    % Nombre de noeuds
Nelt = size(smesh.Elt, 1);        % Nombre d'elements

Ndof = smesh.Ndof_per_node;      % Nombre de ddl par noeud

K = zeros(Nnodes*Ndof, Nnodes*Ndof);

for ee = 1:Nelt
   
   ind = Ndof*repmat(smesh.Elt(ee,2:end),Ndof,1) + repmat((0:Ndof-1)', 1, Ndof) - Ndof + 1;

   ind = ind(:);
   
   K(ind, ind) = K(ind, ind) + ke;  
end

K = sparse(K);