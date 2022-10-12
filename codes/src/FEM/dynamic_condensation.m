function Dc = dynamic_condensation(K, M, frequence, pos_dof_mes, pos_dof_nmes)  

om = 2*pi*frequence;

D = K - om^2*M;

Dmm = D(pos_dof_mes,pos_dof_mes);
Dnmnm = D(pos_dof_nmes,pos_dof_nmes);
Dmnm = D(pos_dof_mes,pos_dof_nmes);

Dnmm = Dmnm.';   
Dc = Dmm - Dmnm*(Dnmnm\Dnmm);