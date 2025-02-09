H = [0, 0, 0];

% Caso 1. Todos los vectores la misma norma
% Prueba definitiva de la ratio
% E_i = [-3, 0.9, 0]

% Caso 2. Todos los vectores la misma norma
% E_i = [-3, 0.9, 1.8]
 
 %Caso 3. Rotar caso 2
 %E_ii = [-3, 0.9, 1.8];
 %ang = 0.37*pi;
 %E_i = Rota_puto_vector(E_ii, ang)
 
 % Caso 4. Otro vector para ver si cambia la ratio
% Pues si que cambia
% E_i = [6, 4, -3]
 
 % Caso 5. Mas sobre la ratio
 % Con 1 vale 1/2, con 2 vale 1 y con 4 vale 2 (la ratio)
 %E_i = [4, 0, 0]
 
 % Caso 6. Todavia mas sobre la ratio
%  E_i = [3, -4, 0]

% Caso 7. Prueba definitiva de la ratio
% La ratio se va a tomar por culo pero si E y H son paralelos, las normas
% de los vectores se mantienen
% E_i = [-3, 0.9, 0];
% H = 0.4 * E_i;

% Caso 8. Vectores perpendiculares
% E_i = [2, 0, 0];
% H = [0, 0.8, 0];
 
% Caso 9. Caso raro para comparar con el mismo rotado
 E_i = [2, -1.3, 3];
 H = [1, 2.5, -0.7];

% Caso 10. Invarianza de la norma total frente a rotaciones
% E_ii = [2, -1.3, 3];
% H_ii = [1, 2.5, -0.7]; 
% ang = 0.21 * pi
% E_i = Rota_puto_vector(E_ii, ang)
% H = Rota_puto_vector(H_ii, ang)
 
% Caso 11. Vectores perpendiculares no sobre el eje
% Los vectores siguen siendo ortogonales
% E_ii = [2, 0, 0];
% H_ii = [0, 0.8, 0]; 
% ang = 0.18*pi;
% E_i = Rota_puto_vector(E_ii, ang)
% H = Rota_puto_vector(H_ii, ang)
 
 
 norm_Ei_H = norm(E_i) + norm(H)
 

ener = 0.5 * norm(E_i)^2 + 0.5 * norm(H)^2
 
 Js_i = Corriente_del_momento(E_i, H)
 
 
 norm_Jx = norm(Js_i(:,1) )
 
 norm_Jy = norm(Js_i(:,2) )
 
 norm_Jz = norm(Js_i(:,3) )
 
 norm_Tot = sqrt( norm_Jx^2 + norm_Jy^2 + norm_Jz^2 )
 
 Ratio = norm_Jx / norm_Ei_H;
  
 % Debe valer 1
 Comprob = 2 * Ratio / norm_Ei_H;
  
Comprob2 = norm_Jx / ener;


% Ver si los vectores son ortogonales
prod_XY = dot(Js_i(:,1), Js_i(:,2)) / norm(Js_i(:,1)) / norm(Js_i(:,2)) 
prod_XZ = dot(Js_i(:,1), Js_i(:,3)) / norm(Js_i(:,1)) / norm(Js_i(:,3)) 
prod_YZ = dot(Js_i(:,2), Js_i(:,3)) / norm(Js_i(:,2)) / norm(Js_i(:,3)) 
 
 