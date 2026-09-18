% La idea de este Matlab es que hay un MCU a velocidad relativista, pero
% luego hay un observador ajeno, que existe en una posición concreta, y
% está mirando a nuestra partícula y quiere ver si la velocidad del objeto
% en MCU evoluciona de forma coherente tras un periodo entero aplicando su
% precesión (se espera que sí) y además comparar dicha precesión  con la
% precesión que ve el observador central para ver si es la misma o no (lo
% cual no lo sé de entrada).

% El observador se llama `U` y la partícula `A`, `C` es el centro del MCU

% CONCLUSION: LA PRECESIÓN MARCADA POR EL SRI CENTRAL DEBE SER APLICADA AL
% RECALCULO DE VELOCIDADES DE CUALQUIER OTRO SRI.
% ESTO ES AL MENOS UN PEQUEÑO EUREKA!

% addpath('C:\Fisica\Scripts matlab\HacerTuTeoria\8_origenes')

R = 1;
% Esta distancia está vista desde `U`.
num_iters = 10000;

v = 0.4;
% Prueba de control:
%v_u  = [0, 0, 0];
% Prueba de interés:
%v_u  = [0, 0.5, 0];

% Prueba máss dura (tambien va bien):
v_u  = [0, -0.2, 0.4];

gamma = fGamma(v);

% Frecuecia angular en el SRI local de la partícula en MCU
w_prim = gamma * v / R;

T_prim = 2 * pi / w_prim;
dt_prim = T_prim / num_iters;

fase_ini = 0;
fase = fase_ini;
dif_fase = 2 * pi / num_iters;

% w_prim^2 o gamma^2 * w^2 pero no ambas cosas!!
a_abs = w_prim^2 * R;

v_ac_ini = - v * [-sin(fase_ini), cos(fase_ini), 0];

% Tiene que estar alineada con `u` para que la transformación de
% velocidades inicial no lo joda todo!
v_au_ini = Vel_Addition_Law(v_u, -v_ac_ini)
v_au = v_au_ini;

iter = 0;
precesion_acum = 0;

% Obtenido de la prueba de control, para debug (obviamnete si v cambia,
% este valor ha de cambiar.
precesion_total_centro = 0.5723;
angulo_chungo = precesion_total_centro / num_iters;

while iter < num_iters
    
    if iter == 5000
        abd = 2131;
    end
    
    a = -a_abs * [cos(fase), sin(fase), 0];
    dv = a * dt_prim;
    
    % Ahora hacer un boost a la velocidad con la que `A` ve a `U`.
    v_au_new = Vel_Addition_Law(v_au, dv);
    
    % Probaturas chungas: ¡RESULTA QUE DAN LA VICTORIA!
    mat_rot = MatrizRotacion([0,0,-1], angulo_chungo);
    
    %[dir_rot, angulo] = Angulo_Wigner(v_au_new, dv);
    % El angulo tiene que ser negativo (o la matriz transpuesta)
    %mat_rot = MatrizRotacion(dir_rot, angulo);
    
    % 18.09. Veamos si esto se puede hacer con la matriz que va mejor
    % No usamos más el angulo chungo, podemos usar -v_ca
    v_ca = v * [-sin(fase), cos(fase), 0];
    RW = RotacionWigner(-v_ca, dv);
    RW_3 = RW(2:4, 2:4);
    
    v_au_new = v_au_new * mat_rot';
    %v_au_new = v_au_new * RW_3';
    
    v_au = v_au_new;
    
    %precesion_acum = precesion_acum + angulo;
    precesion_acum = precesion_acum + angulo_chungo;    
    
    fase = fase + dif_fase;
    iter = iter + 1; 
end

fase_final = fase - 2 * pi
precesion_total = precesion_acum
v_au_final  = v_au
