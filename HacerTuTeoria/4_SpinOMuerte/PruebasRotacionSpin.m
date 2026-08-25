% Script para ver cómo cambia el spin con la rotación, para ver qué es eso
% de que el spin ha de cambiar 720º para volver a su posición inicial.

% p=0, sistema en reposo!!!
 p0 = [0, 0, 0];


MsZ = MatrizMs(3,1);         % Y
ph_ini = [1, 0];     %  Z
Psi_ini = DiracSpinorPlainWave(p0, ph_ini)
Psi = Psi_ini;

dfase = pi / 1000;
n = 0;
%while n < 500      % [0.5, 0.5, 0, 0]/raiz(2);
while n < 1000     % [0, 1, 0, 0];
%while n < 1500     % [-0.5, 0.5, 0, 0]/raiz(2);
%while n < 2000     % [-1, 0, 0, 0];
%while n < 2500     % [-0.5, -0.5, 0, 0]/raiz(2);
%while n < 3000     % [0, -1, 0, 0];
%while n < 3500     % [0.5, -0.5, 0, 0]/raiz(2);
%while n < 4000      % [1, 0, 0, 0]/raiz(2);

    dPsi = dfase * SpinorRotation(MsZ, Psi)';
    Psi = Psi + dPsi;
    
    Psi = Psi / norm(Psi);
    
    n = n + 1;
end

fase_fin = dfase * n
Psi = Psi


% Otra prueba: ver cómo cambia la fase global rotando en el eje del spin

% MsZ = MatrizMs(1,2);         % Y
% ph_ini = [1, 0];     %  Z
% Psi_ini = DiracSpinorPlainWave(p0, ph_ini)
% Psi = Psi_ini;
% 
% dfase = pi / 1000;
% n = 0;

%while n < 500               % [1 + i, 0, 0, 0]/raiz(2)
%while n < 1000              % [i, 0, 0, 0]
%while n < 1500              % [-1 + i, 0, 0, 0]/raiz(2)
%while n < 2000              % [-1, 0, 0, 0]
%while n < 2500              % [-1 - i, 0, 0, 0]/raiz(2)
%while n < 3000              % [-i, 0, 0, 0]
%while n < 3500              % [1 - i, 0, 0, 0]/raiz(2)
%while n < 4000              % [1, 0, 0, 0]
% 
%     dPsi = dfase * SpinorRotation(MsZ, Psi)';
%     Psi = Psi + dPsi;
%     
%     Psi = Psi / norm(Psi);
%     
%     n = n + 1;
% end
% 
% fase_fin = dfase * n
% Psi = Psi
