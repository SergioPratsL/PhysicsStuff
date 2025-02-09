
function v = ObtenVelocidadDeMatrizBoost(matrizBoost)
% matrizBoost es una matriz 4x4 que representa un boost
% Esta función no comprueba que matrizBoost sea un boost, si le envías una
% patata, simplemente te devolverá otra

v = matrizBoost(1,2:4) / matrizBoost(1,1);