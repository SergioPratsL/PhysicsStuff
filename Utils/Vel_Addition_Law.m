function v_otro_en_new_SRI = Vel_Addition_Law( v_otro, v_new_SRI );
% Ambas velocidades normalizadas a la velocidad de la luz.
% Nombres suficientemente intuitivos

% 15/05/2018. Esto es una puta mierda que solo sirve para direcciones
% paralelas!
%v_otro_en_new_SRI = (v_otro - v_new_SRI) / (1 - dot(v_otro, v_new_SRI));

vector = [1, v_otro];
vector_boost = Boost(vector, v_new_SRI);

v_otro_en_new_SRI = vector_boost(2:4) / vector_boost(1);
