clear;
% Como era de esperar, esto no resolvió nada.
% 
% Sx = MatrizSpin_4_4(1);
% 
% phi = sym('phi', [4 1]);
% phi_dag = conj(phi.');
% 
% rho = phi_dag * phi;
% phi_rot = phi_dag * Sx * phi;
% 
% terminoChungo = (phi_rot)/rho * phi;
% 
% pretty(terminoChungo)



Sx = MatrizSpin_4_4(1);
[gt, gx, gy, gz] = MatricesGamma();
[at, ax, ay, az] = MatricesAlfa();
g5 = MatrizGamma(5);

%Pruebas:
v = [0, 0, 0];
%spinor_base = [1, 0]; 
spinor_base = [1, 1]; 
%spinor_base = [cos(pi/8), sin(pi/8)]    % 45!


%v = [0.4, 0, 0];
v = [0, 0.4, 0];

vector_original = SpinorToVector(spinor_base)

spinor_base = spinor_base / norm(spinor_base);
p = fGamma(v) * v;

phi = DiracSpinorPlainWave(p, spinor_base).'

rotation = -i * Sx * phi
proyeccionRotacion = phi' * rotation 

% Parece que pretender sumar el vector de rotación es un error, porque va
% mal, pero si lo hago directamente con el spin va bien.
vector_rotacion = SpinorToVector(rotation(1:2).')

vector_rotado = SpinorToVector(phi(1:2).' + 0.01 * rotation(1:2).')



%ojalaQueRote = Sx * phi - Sx * gx * phi
%ojalaQueRote = Sx * phi - Sx * ax * phi
%ojalaQueRote = Sx * phi - g5 * ax * phi   % Se lo carga todo!
ojalaQueRote = Sx * phi - Sx * gx * phi;

comprobacion = phi' * ojalaQueRote;

