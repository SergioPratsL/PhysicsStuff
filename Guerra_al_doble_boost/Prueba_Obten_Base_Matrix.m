clear;

eje = [0,0,1];

B = ObtenBaseMatrix(eje)

result = (B * eje')'

ang = pi / 12;

matrix_base = [1, 0, 0; 0, cos(ang), -sin(ang); 0, sin(ang), cos(ang)];

matrix = B' * matrix_base * B



