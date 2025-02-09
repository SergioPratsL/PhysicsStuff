function E_o = Rota_puto_vector(E_i, ang)

E_o = E_i;

E_o(1) = cos(ang) * E_i(1) - sin(ang) * E_i(2);
E_o(2) = cos(ang) * E_i(2) + sin(ang) * E_i(1);

