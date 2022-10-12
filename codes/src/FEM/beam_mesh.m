function bmesh = beam_mesh(xmin, L, Nelt)

Nnodes = Nelt + 1;

xmax = L + xmin;

ID = 1:Nnodes;
x = xmin: L/Nelt : xmax;

bmesh.Nodes = [ID(:), x(:)];

bmesh.Elt = zeros(Nelt,3);

for ee = 1:Nelt
   bmesh.Elt(ee,:) = [ee, ee, ee+1];
end

bmesh.Ndof_per_node = 2;
bmesh.elem_size = L/Nelt;