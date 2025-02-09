
R_vect = [1,0,0];
v_vect = [0.4, 0, 0];
despl = [1, 3, 1];
ac_emis = [0,1,0];
dv_recept = despl;


R_vect = [1,-1.5,0];
v_vect = [0.4, 0, 0];
despl = [1, 3, 1];
ac_emis = [0,1,0.8];
dv_recept = despl;


R_vect = [1,-1.5,4];
v_vect = [0.2, 0.18, 0];
despl = [1, 3, 1];
ac_emis = [-2,1,0.8];
dv_recept = despl;



Dif_Pot_B_EM = EfectosCamposRadiado(R_vect, v_vect, despl, ac_emis, dv_recept)