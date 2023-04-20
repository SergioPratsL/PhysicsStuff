function vect_rot = ThomasPrecession( Sigma, dv, v, v_new );
% Sigma está asociado a v, se puede pasar para ahorrar operaciones.
% Recordar, dv debe ser variación de velocidad desde el SRI del lab.
% 09.01.2023 El punto v_new lo usé por algún motivo, para no perder el
% cambio a v_new, si no hace falta pasar lo mismo que en v.

if Sigma == 0
    Sigma = fSigma(v);
end

[eje_rot, angulo] = ThomasPrecessionAng( Sigma, dv, v );

rotMatrix = RotationMatrix(angulo);

% rotMatrix se traspone para rotar en sentido contrario.
vect_rot = RotacionGeneral(v_new, eje_rot, rotMatrix');

