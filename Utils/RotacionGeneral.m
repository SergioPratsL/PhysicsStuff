function vector_rotado = RotacionGeneral( vector_a_rotar, eje_rot, rotMatrix );
% Rota un vector espacial sobre cualquier eje
% Recibe la matriz de rotación en vez del ángulo a rotar parar ahorrar
% calculos

[e1, e2, e3] = ObtenBase(eje_rot);
vector_new_base = ConvierteNuevaBase(vector_a_rotar, e1, e2, e3);
vector_new_base_rotado = vector_new_base * rotMatrix;
vector_rotado= ConvierteBaseOriginal(vector_new_base_rotado, e1, e2, e3);