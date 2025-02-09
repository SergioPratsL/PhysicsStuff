function matrix = RotationMatrix( ang )
% Rotación sobre el eje X, angulo en radianes!

matrix = [1, 0, 0; 0, cos(ang), -sin(ang); 0, sin(ang), cos(ang)];

end