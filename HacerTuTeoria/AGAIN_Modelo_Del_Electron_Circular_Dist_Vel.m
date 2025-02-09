% Este modelo sólo tiene dos elementos pero la gracia es que uno va a
% diferente velocidad

me = 9.1093 * 10^-31;
c = 299792458;
h = 6.62607 * 10^-34;
h_bar = h /(2*pi);
qe = 1.6022 / 10^(19);
perm_elec = 8.854 *  10^-12;
E_elec = me * c^2;

r_compton = h/(c*me);

% Con este radio parece que va, pero es muy raro... es un tinglado.
r = r_compton;

% factor multiplicador
%a = 1;
% a = 0.11727;  % Usando el momento alternativo.
a = 0.08775;
r = r * a;

cof2 = 1 - 1/(a*4*pi);
cof1 = 1;
cof0 = -1/(a*4*pi);



p = [cof2, cof1, cof0];
sols = roots(p);

% Espero que cambie...
v2 = 2;

if norm(imag(sols(1))) < 10^-8 && norm(sols(1)) <= 1
    v2 = sols(1)
end

if norm(imag(sols(2))) < 10^-8 && norm(sols(2)) <= 1
    v2 = sols(2)
end

if( abs(v2) > 1)
    disp('velocidad incoherente!')
    v2 = v2
    return
end

% Recordar: m1 y m2 son coeficientes, deben multiplicarse por me.
m1 = 1 * v2^2 / (1+v2^2)
m2 = 1 - m1

if( m1 < 0 || m1 > 1)
    disp('masas incoherentes')
    return
end

q1q2 = - c^2 * r /qe^2  * (4*pi*perm_elec*me*m1)


cof2 = 1;
cof1 = 1;
cof0 = q1q2;

p2 = [cof2, cof1, cof0];
sols2 = roots(p2)

q1 = sols2(1);
q2 = sols2(2);

mom_mag = qe*c*(q1 +v2*q2)*r
mom_mag_altern = qe*c*(q2 +v2*q1)*r;

mom_mag_esperado = -2 * (qe/me)*(h_bar/2)

L = r*c*(m1*me + v2*m2*me)

L_rel_hbar = L / h_bar


% último check:
F_1 = (q1*q2*qe^2)/(4*pi*perm_elec)/r^2;
% El signo de la aceleración es bueno porque el signo de q1 y q2 es
% diferente
acel_1_abs = abs(F_1 /(me*m1))

acel_1_necesaria = c^2 / r



r_basico = h_bar / (2*c*me)

r_basico_div_r_comton = r_basico / r_compton