
eje_rot = [0,0,1];
w = pi/6;

RotMatrix = RotationMatrixGeneral(eje_rot, w)

% Siendo en el mismo eje, las rotaciones se acumulan
%RotMatrix2 = RotMatrix * RotMatrix;

% Equivale a un giro de 90º entre la X y la Z.
MiNuevaBase = [0,0,-1;0,1,0;1,0,0]

% Ya que ahora la X está en la tercera columna, deberían rotar Y y Z (la Z
% en lugar de la X).
RotMatrixNuevaBase = MiNuevaBase' * RotMatrix * MiNuevaBase

determinante = det(RotMatrixNuevaBase)

% Todo parece bien, se obtiene la rotación en la otra base

% No esto no está bien!
%RotarEnBaseOriginal = RotMatrixNuevaBase * MiNuevaBase'

% Tiene que hacerse todo esto para volver a la base original!
RotarEnBaseOriginal = MiNuevaBase * RotMatrixNuevaBase * MiNuevaBase'


% Si sólo haces esto: MiNuevaBase' * RotMatrix
% lo que obtienes es una matriz ortonormal de determinante 1 pero parece
% que la diagonal se ha girado y coge el 31-22-13 (o 13-22-31 según se
% mire)

% Ahora giraré 45º sobre el plano YZ (en origen)
eje_rot2 = sqrt(1/2)*[1,0,0];
w2 = pi/4;

RotMatrix2 = RotationMatrixGeneral(eje_rot2, w2)

RotMatrixNuevaBase2 = RotMatrixNuevaBase' * RotMatrix2 * RotMatrixNuevaBase;

X = RotMatrix * RotMatrix2
Y = RotMatrix2 * RotMatrix
