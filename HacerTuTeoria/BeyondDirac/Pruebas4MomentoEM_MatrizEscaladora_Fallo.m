
clear;

%Prueba 1
%E = [3, 0, 0];
%B = [0, 1, 0];


%Prueba Final
E = [11, 13, 17];
B = [3, 5, 7];


u = norm(E)^2 + norm(B)^2
p = cross(E, B)

F = TensorEM(E, B);

[gt, gx, gy, gz] = MatricesGamma();

F_trasp = F'
cosa = F_trasp * gt


M_izq = F' * gt;
M_der_x = gx * F;
M_der_y = gy * F;
M_der_z = gz * F;
M_der_t = gt * F;
%M_izq = F';
%M_der_t = -F / 2;

px_calc = -1i * CompactaMatrices(M_izq, M_der_x) / 2
py_calc = -1i * CompactaMatrices(M_izq, M_der_y) / 2
pz_calc = -1i * CompactaMatrices(M_izq, M_der_z) / 2
pt_calc = -CompactaMatrices(M_izq, M_der_t) / 2


% M_izq_x = F' * gx;
% M_izq_y = F' * gy;
% M_izq_z = F' * gz;
% M_izq_t = F' * gt';
% 
% px_calc = CompactaMatrices(M_izq_x', M_der_x)
% py_calc = CompactaMatrices(M_izq_y', M_der_y)
% pz_calc = CompactaMatrices(M_izq_z', M_der_z)
% pt_calc = CompactaMatrices(M_izq_t', M_der_t)
% Obviamente todas dan lo mismo salvo por el signo...


function val = CompactaMatrices(M1, M2)
    val = 0;

    m = 1;
    while m < 5
        n = 1;
        while n < 5
            val = val + M1(m,n) * M2(m,n);
            n = n + 1;
        end

        m = m + 1;
    end
end