% La mía de lo mismo que la de gemini, pues me quedo con la mía.
% 12.09.2026. Pues no iba tan bien para precesiones parece que sí, pero
% para angulos grandes no!
% 17.09.2026. Pues la otra funcion tampoco va bien del todo así que me
% quedo con esta, pero esto sólo sirve para calcular precesiones cuando
% tanto v1 como v2 son grandes, usar RotacionWigner.

function [dir_rot, angulo] = Angulo_Wigner( v1, v2 )

gam1 = fGamma(v1);

gam2 = fGamma(v2);

dir_rot = cross(v1, v2);
dir_rot = dir_rot / norm(dir_rot);

tanRho = gam1 * gam2 * norm(cross(v1, v2)) / (gam1 + gam2);

angulo = atan(tanRho);

% function [axis, angle] = Angulo_Wigner(v1, v2)
%     Angulo_Wigner Computes Wigner rotation angle (rad) and axis
%     for consecutive boosts v1 followed by v2 (c=1).
%     
%     Inputs:
%       v1, v2 - 3-element velocity vectors (normalized to c=1)
%     
%     v1 = v1(:); 
%     v2 = v2(:);
%     
%     b1_sq = dot(v1, v1);
%     b2_sq = dot(v2, v2);
%     
%     if b1_sq >= 1 || b2_sq >= 1
%         error('Velocity magnitudes must be strictly less than 1.');
%     end
% 
%     Cross product determines non-collinearity
%     cp = cross(v1, v2);
%     cp_norm = norm(cp);
%     
%     if cp_norm < 1e-12 || b1_sq < 1e-12 || b2_sq < 1e-12
%         axis = [0; 0; 1];
%         angle = 0;
%         return;
%     end
% 
%     Axis is always perpendicular to the plane formed by v1 and v2
%     axis = cp / cp_norm;
% 
%     Gamma factors
%     g1 = 1 / sqrt(1 - b1_sq);
%     g2 = 1 / sqrt(1 - b2_sq);
% 
%     Exact analytical formula derived from trace of Wigner rotation matrix
%     tan(theta / 2) = |v1 x v2| * g1 * g2 / (1 + g1 + g2 + g1 * g2 * (1 + dot(v1, v2)))
%     numerator = cp_norm * g1 * g2;
%     denominator = 1 + g1 + g2 + g1 * g2 * (1 + dot(v1, v2));
%     
%     half_angle = atan(numerator / denominator);
%     angle = 2 * half_angle;
% end