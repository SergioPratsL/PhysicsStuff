
function [Var1, Var2, Var3] = PruebaFuego1( R, v1, v2 )

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Esta es la primera prueba de fuego para el potencial absoluto

% Var1 es el campo electrico perpendicular a al plano Pl
% formado por R y (v1-v2) y al menos para v1 = 0 o v2 = 0 debe valer 0

% Var2 es la variacion del potencial por unidad de movimiento en la
% direccion de la componente de (v1-v2) (La X en mis notas).

% Var3 es la variacion del potencial por unidad de movimiento en la
% direccion de R medida segun la componente del plano Pl ortogonal a
% (v1-v2)

% Logicamente, Var2  y Var3 deben valer lo mismo

% Que se cumplan estas premisas es la primera prueba de fuego

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% ... Lo que se evalua son los campos vistos por "2"

% R es la distancia medida de "2" a "1" (en el diario lo hice asi...)
% v1 es la velocidad con la que se mueve "1", primero probar que sea 0
% v2 es la velocidad con la que se mueve "2" (que no sea 0, maldito cobarde!!)

% Lo defino desde el punto de vista de  2.
dif_v = v2 - v1; 

% Ojo, el "Dv" donde la D deberia ser un triangulo, puede ser positivo o
% negativo pero esta direccionalidad la contiene "dif_v_norm"
Dv = norm(dif_v);

% Prueba trivial!
if Dv == 0
    return
end

dif_v_norm = dif_v / Dv;

R_norm = R / norm(R);

% Vector que es ortogonal al plano Pl
vector_ortogonal = cross( dif_v_norm, R_norm );

norm_v_ort = norm(vector_ortogonal);

if  norm_v_ort > 0
    
    vector_ortonormal = vector_ortogonal / norm(vector_ortogonal);

% En este caso la distancia y la diferencia de lecturas estaban alineadas!!
% elegir el vector que me salga del rabo!
else

    if dif_v_norm(2) ~= 0  ||  dif_v_norm(3) ~= 0
        vect_test = [1, 0, 0];
    else
        vect_test = [0, 1, 0];
    end
    
    vect_ortogonal = cross( dif_v_norm, vect_test );
    
    vector_ortonormal = vector_ortogonal / norm(vector_ortogonal);    
    
end

% Obtenemos el segundo vector del plano Pl
vector_planar = cross( dif_v_norm, vector_ortonormal );
vector_planar = vector_planar / norm(vector_planar);

% Los famosos valores de las formulas, Rx y Ry
% Entiendase que "x" es la direccion de "v2-v1"
% mientras que "y" es la direccion perpendicular 
% dentro del plano que forma "v2-v1" y R.
Rx = dot( R, dif_v_norm);

Ry = dot( R, vector_planar);


% Obtener las distancias retardadas y los cambios EM de "1" sobre "2" y de
% "2" sobre "1"

Rr12 = DRetardada(-R, v1);
Rr21 = DRetardada(R, v2);

% Paso de las unidades ya que me la liaran parda, prefiero trabajar con
% numeros proximos a 1...
[E12, B12] = CampoInducido_sin_unidades(Rr12, v1);
[E21, B21] = CampoInducido_sin_unidades(Rr21, v2);

% Rotar todos los vectores
E12_rot_x = dot(E12, dif_v_norm);
E12_rot_y = dot(E12, vector_planar);
E12_rot_z = dot(E12, vector_ortonormal);
E12_rot = [E12_rot_x, E12_rot_y, E12_rot_z];

B12_rot_x = dot(B12, dif_v_norm);
B12_rot_y = dot(B12, vector_planar);
B12_rot_z = dot(B12, vector_ortonormal);
B12_rot = [B12_rot_x, B12_rot_y, B12_rot_z];

E21_rot_x = dot(E21, dif_v_norm);
E21_rot_y = dot(E21, vector_planar);
E21_rot_z = dot(E21, vector_ortonormal);
E21_rot = [E21_rot_x, E21_rot_y, E21_rot_z];

B21_rot_x = dot(B21, dif_v_norm);
B21_rot_y = dot(B21, vector_planar);
B21_rot_z = dot(B21, vector_ortonormal);
B21_rot = [B21_rot_x, B21_rot_y, B21_rot_z];

% Tb hay que rotar las velocidades o sera un desastre :P
v1_rot_x = dot(v1, dif_v_norm);
v1_rot_y = dot(v1, vector_planar);
v1_rot_z = dot(v1, vector_ortonormal);
v1_rot = [v1_rot_x, v1_rot_y, v1_rot_z];

v2_rot_x = dot(v2, dif_v_norm);
v2_rot_y = dot(v2, vector_planar);
v2_rot_z = dot(v2, vector_ortonormal);
v2_rot = [v2_rot_x, v2_rot_y, v2_rot_z];


Vect_Var1 = cross( v2_rot, B12_rot ) + cross( v1_rot, B21_rot);

%Var1 = dot( Vect_Var1, vector_ortonormal );
% En la nueva base Var1 es la tercera componente de Vect_Var1
Var1 = Vect_Var1(3);

Prod_E_v = dot(E12_rot, v2_rot) + dot(E21_rot, v1_rot);

Vect_aux_2_3 = E21_rot + cross( v2_rot, B12_rot ) + cross( v1_rot, B21_rot);

Escalar_aux_2 = (1/Dv) * (dot(E12_rot, v2_rot) + dot(E21_rot, v1_rot) );

Var2 = Vect_aux_2_3(1) + Escalar_aux_2;

Escalar_aux_3 = Escalar_aux_2 * (Rx/Ry);

Var3 = - (Ry / norm(R)) * ( Vect_aux_2_3(2) - Escalar_aux_3);






