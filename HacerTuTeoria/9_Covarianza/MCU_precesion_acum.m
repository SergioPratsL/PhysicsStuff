% Esta version de MCU_visto_por_otro.m es para ver si la precesión crece
% mucho cuando la velocidad del MCU crece mucho

R = 1;
% Esta distancia está vista desde `U`.
num_iters = 40000;

v = 0.95;

gamma = fGamma(v);

% Frecuecia angular en el SRI local de la partícula en MCU
w_prim = gamma * v / R;

T_prim = 2 * pi / w_prim;
dt_prim = T_prim / num_iters;

fase_ini = 0;
fase = fase_ini;
dif_fase = 2 * pi / num_iters;

a_abs = w_prim^2 * R;

v_ac_ini = - v * [-sin(fase_ini), cos(fase_ini), 0];

iter = 0;
precesion_acum = 0;
v_ac = v_ac_ini;

while iter < num_iters
    
    if iter == 5000
        abd = 2131;
    end
    
    a = -a_abs * [cos(fase), sin(fase), 0];
    dv = a * dt_prim;
    
    v_ac_new = Vel_Addition_Law(v_ac, dv);
    
    %[dir_rot, angulo] = Angulo_Wigner(v_ac, dv);
    %% El angulo tiene que ser negativo (o la matriz transpuesta)
    %mat_rot = MatrizRotacion(dir_rot, angulo);    
    
    RW = RotacionWigner(v_ac, dv);
    RW_3 = RW(2:4, 2:4);
    
    [eje_rot, angulo] = ExtraeAnguloDeMatriz(RW_3);

    v_ac_new = v_ac_new * RW_3';
    
    v_ac = v_ac_new;
    
    precesion_acum = precesion_acum + angulo;
    
    fase = fase + dif_fase;
    iter = iter + 1; 
end

fase_final = fase - 2 * pi
precesion_total = precesion_acum * 180 / pi
v_ac_final  = v_ac
