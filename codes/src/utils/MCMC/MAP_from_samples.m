function Fmap = MAP_from_samples(F)

Np = size(F, 1);

Fmap = zeros(Np, 1);

for ee = 1:Np
   [pdf, x] = kernelDensity(F(ee,:));
   
   Fmap(ee) = x(pdf == max(pdf));
end