syms a

x_2_0 = 3 * cos(a)^2-1;

dx_2_0 = diff(x_2_0,a);
ddx_2_0 = diff(dx_2_0,a);

% dRho^2(Ond)*Ond*sin(Rho)
x_2_0_int = ddx_2_0 * sin(a) * x_2_0;

y_2_0 = int(x_2_0_int, a);

% x7 = 6*(sin(a)^2-cos(a)^2)*(3*cos(a)^2-1)*sin(a)
% y7 = int(x7, a)

% El resultado es 6*cos(a) - 10*cos(a)^3 + (36*cos(a)^5)/5
% Al integrar este valor por los dos hemisferios y por el factor 2pi de
% longitud de cada "paralelo" y multiplicar por (5/16pi), el resultado es
% 4, o sea, 2^2.


x_1_0 = cos(a);
dx_1_0 = diff(x_1_0, a);
ddx_1_0 = diff(dx_1_0, a);

x_1_0_int = ddx_1_0 * sin(a) * x_1_0;

y_1_0 = int(x_1_0_int, a);

% Ahora el verdadero L^2

% sin(a) se simplifica con 1/sin(a) 
ddx_1_0_esp = diff( (sin(a)*dx_1_0), a);
x_1_0_int_esp = ddx_1_0_esp * x_1_0;

y_1_0_esp = int(x_1_0_int_esp, a)


ddx_2_0_esp = diff( (sin(a)*dx_2_0), a);
x_2_0_int_esp = ddx_2_0_esp * x_2_0;

y_2_0_esp = int(x_2_0_int_esp, a)


