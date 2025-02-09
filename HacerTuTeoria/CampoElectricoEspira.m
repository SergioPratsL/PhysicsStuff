clear;

% La carga q=1, eso no se discute aquí.

% Quiero < 1 para dar vidilla a la aceleración
radio = 1;

% Velocidad a la que las cargas se mueven (en unidades naturales).
v_mod = 0.1;

num_iteraciones = 360;
dif_fase = 2 * pi / num_iteraciones;

iter = 0;


%r_eval = [3,0,0]
% E =   1.0e-14 *   [-0.0333    0.1164         0]
% B =    0         0   -0.7609

r_eval = [1.41,0,2]
% E = 1.0e-14 * [-0.0923    0.0917    0.1443]
% B =  1.4320   -0.0000    1.2936


% Confirmado, E=0 en la espira de corriente

v_0 = [0,0,0];

E = 0;
B = 0;
E_prot = 0;

while iter < num_iteraciones
    fase = iter * dif_fase;
   
    r0 = radio * [cos(fase), sin(fase), 0];
    
    v = v_mod * [-sin(fase), cos(fase), 0];
    
    a = - (v_mod^2 / radio) * [cos(fase), sin(fase), 0];
    %a = v_0;
    
    r_ret = r_eval - r0;
    
    if iter == 90
        iter = iter;
    end
    
    Dop_Esp = 1 - dot(v,r_ret)/norm(r_ret);
    
    [e_ind, b_ind] = CampoInducido_sin_unidades( r_ret, v);
    [e_rad, b_rad] = CampoRadiadoCargaAcelerada( r_ret, v, a);
    
    % Descontar la parte estática
    [e_prot, b_prot] = CampoInducido_sin_unidades( r_ret, v_0);
    
    
    E = E + (e_ind + e_rad) * Dop_Esp - e_prot;
    B = B + (b_ind + b_rad) * Dop_Esp;
    %E_prot = E_prot + e_prot;
      
    iter = iter + 1;
end

E = E
B = B
S = cross(E, B)

%E_prot = E_prot


