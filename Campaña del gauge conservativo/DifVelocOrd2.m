function [dv_dvrdve, dv_dvedvr] = DifVelocOrd2( v_vect, dv_e, dv_r );
% Package: aceleracion dual

% nomenclatura: la ultima derivada es la que esta mas cerca de a velocidad
%dv_dvrdve es  d(dv/dve)/dvr
%dv_dvedvr es  d(dv/dvr)/dve

% NOTA: EL PROBLEMA ES QUE LAS COMPONENTES PERPENDICULAR Y PARALELA TAMBIEN
% CAMBIAN

    
v = norm(v_vect);
v_norm = v_vect / v;

Sigma = 1 / sqrt(1-v^2);

dv_e_paral = v_vect * dot(dv_e, v_norm);
dv_e_perp = dv_e - dv_e_paral;


dv_r_paral = v_vect * dot(dv_r, v_norm);
dv_r_perp = dv_r - dv_r_paral;

dv_dvrdve = dot(v_vect, dv_r) / Sigma * dv_e_perp + 2 * dot(v_vect, dv_r) * dv_e_paral;

dv_dvedvr = dot(v_vect, dv_e) / Sigma * dv_r_perp + 2 * dot(v_vect, dv_e) * dv_r_paral;