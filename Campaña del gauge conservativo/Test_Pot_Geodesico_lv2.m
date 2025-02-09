% Package: prueba geodesica nivel 2

% Esta prueba busca, para salvar los trastos, justificar que las
% aceleraciones que se hacen las particulas B y C entre si no modificaran en
% conjutno el potencial que ve A de ambas (la suma de ambos potenciales)

% Mi intencion es plantarme aqui si sale esta prueba y pedir hacer una
% tesis o yo que se que...

clear


% Prueba 1
R_B = [0.5,1,0];
v_B = [0.3, 0, 0];

R_C = [0.5,-1,0];
v_C = [-0.3, 0, 0];




[Dif_Pot_B_B, Dif_Pot_B_A] = VarPotGeodesico_lv2(R_B, v_B, R_C, v_C)

[Dif_Pot_C_B, Dif_Pot_C_A] = VarPotGeodesico_lv2(R_C, v_C, R_B, v_B)


Dif_Tot = Dif_Pot_B_B + Dif_Pot_B_A + Dif_Pot_C_B + Dif_Pot_C_B


Dif_aux_1 = Dif_Pot_B_B + Dif_Pot_C_B

Dif_aux_2 = Dif_Pot_B_A + Dif_Pot_C_A

