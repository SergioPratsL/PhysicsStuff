% EjeRot no es más que un comentario porque son los parámetros de MatrizMs
% los qe determinan cual es el eje de rotación!

% KK
% Ms1 = MatrizMs(2,3);        % X
% p1 = [0, 0, 0.5];
% ph1 = [1, -1] / sqrt(2);     % -X
% 
% Psi = DiracSpinorPlainWave(p1, ph1)
% dPsi = SpinorRotation(Ms1, Psi)


% p=0, sistema en reposo!!!
 p0 = [0, 0, 0];

 
% Ms2 = MatrizMs(2,3);          % X
% ph2 = [1, -1] / sqrt(2);      % -X
% 
% Psi = DiracSpinorPlainWave(p0, ph2)
% dPsi = SpinorRotation(Ms2, Psi)'
% Ratios: i*[-1/2, -1/2, N/A, N/A]
% Básicamente cambia la fase, en el mismo sentido en las dos componentes.



% Ms3 = MatrizMs(2,3);        % X
% ph3 = [1, +1i] / sqrt(2);   % +Y
% 
% Psi = DiracSpinorPlainWave(p0, ph3)
% dPsi = SpinorRotation(Ms3, Psi)'
% Ratios: i*[-1/2, 1/2, N/A, N/A]
% Se mueve hacia la Z, lógico para un Spin en Y rotando sobre X :).



% Ms4 = MatrizMs(2,3);        % X
% ph4 = [1, 0];               % +Z
% 
% Psi = DiracSpinorPlainWave(p0, ph4)
% dPsi = SpinorRotation(Ms4, Psi)'
% No hay ration dPSi parece que se orienta hacia la segunda componente, y
% en imaginario, pero eso da igual porque está vacía, es un giro del spin
% hacia Y. Si giras sobre el eje X una cosa en Z, tiene que ir hacia Y.


Ms5 = MatrizMs(1,2);          % Z
ph5 = [1, 0];                 % Z

Psi = DiracSpinorPlainWave(p0, ph5)
dPsi = SpinorRotation(Ms5, Psi)'
% Ratios: i*[1/2, N/A, N/A, N/A], Cambia la fase :)



% Ms6 = MatrizMs(1,2);          % Z
% ph6 = [1, -1i] / sqrt(2);     % -Y
% 
% Psi = DiracSpinorPlainWave(p0, ph6)
% dPsi = SpinorRotation(Ms6, Psi)'
%Ratios: i*[1/2, 1/2, N/A, N/A]  --> Gira hacia -X (lo se por elSpinorToVector).


% Otras pruebas que poco tienen que ver pero que quiero poner aquí:

% ph = [1, 0];   
% Sigma_x = MatrizPauli(1);
% Sigma_y = MatrizPauli(2);
% Sigma_z = MatrizPauli(3);
% 
% val_x = ph * Sigma_x * ph'
% val_y = ph * Sigma_y * ph'
% val_z = ph * Sigma_z * ph'


% sx = MatrizPauli(1);
% sy = MatrizPauli(2);
% sz = MatrizPauli(3);
% 
% prueba1 = sy - i*sx*sz
% prueba2 = sx + i*sy*sz



