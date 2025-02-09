% Este script quiere ver si la precesión de Thomas se boostea igual que el
% momento angular, lo cual reforzaría la idea que la precesión de Thomas
% genera momento angular.
% Trabajará sobre MCUs, porque el momento angular de un punto no es de fiar

clear;

r_lab = 1;
v_lab = 0.6;
w_lab = v_lab / r_lab;
a_lab = v_lab^2 / r_lab;

dir_mom_lab = [0,0,1];
dir_mom_lab = dir_mom_lab / norm(dir_mom_lab);

v_otro_SRI = [0,0,0.6];
%v_otro_SRI = [0,0,0];   % conseguí que el while cuadrara para este caso.

% Puede demostrarse que al integrar sobre la anill el valor de N es 0.
N_lab = [0,0,0];
% Pero quizá haya que ir punto a punto (pues no ayuda mucho que digamos)
%%%N_lab = [0, r_lab,0];

% Cojo un punto cualquiera del MCU, digamos, la posición +Y.
x_MCU_lab = r_lab * [0, 1, 0];
v_MCU_lab = v_lab * [1, 0, 0];
a_MCU_lab = a_lab * [0, -1, 0];

Sigma_lab = fGamma(v_MCU_lab);
[dir, w] = ThomasPrecessionAng(Sigma_lab, a_MCU_lab, v_MCU_lab);

MatrixRotacionDesdeLab = RotationMatrixGeneral(dir, w)

wT_lab = dir * w;

lT_lab = wT_lab * r_lab * Sigma_lab;

% 4 aceleración en el laboratorio
a_u_lab = [0, a_MCU_lab];

% Ahora hay que trasformar este momento al nuevo SRI, se transforma igual
% que el boost EM.
[N_otro_SRI, lT_otro_SRI] = Boost_EM(N_lab, lT_lab, v_otro_SRI);
lT_otro_SRI = lT_otro_SRI




% El otro valor es calcular la precesión de Thomas en el otro SRI a ver si
% da lo mismo
v_MCU_otro = Vel_Addition_Law(v_MCU_lab, v_otro_SRI);
Sigma_otro = fGamma(v_MCU_otro);

%a_MCU_u_otro = Boost(a_u_lab, v_MCU_otro);
% La aceleración es menor por la dilatación del tiempo propio!
%Sigma_u = fGamma(v_otro_SRI);
%a_MCU_otro = a_MCU_u_otro(2:4) / Sigma_u;


% Tendré que averiguar los valores correctos a base de prueba y error.
%rotMatrix = GetThomasRotMatrix(v_MCU_lab, v_otro_SRI );
%a_MCU_otro_rot = a_MCU_otro * rotMatrix;
% Esto no tiene que rotar si el que manda es el SRI remoto!
%%v_MCU_otro_rot = v_MCU_otro * rotMatrix;  


%[dir, w] = ThomasPrecessionAng(Sigma_otro, a_MCU_otro_rot, v_MCU_otro_rot);
%[dir, w] = ThomasPrecessionAng(Sigma_otro, a_MCU_otro, v_MCU_otro);
%wT_otro = dir * w;
% r_lab no cambia en el ejemplo más fácil
%lT_otro_Thomas = wT_otro * r_lab
%lT_otro_Thomas_rot = lT_otro_Thomas * rotMatrix

% ESTO NO VA.

% Voy a intentar una cosa más difícil.
MatrizRotacionTotal = eye(3);

NumSegmentos = 1000;

n = 0;
w_otro = w_lab / fGamma(v_otro_SRI);
wT_otro = 0;

% Empezamos por el 0 pero lo ignoramos deliveradamente.
while n <= NumSegmentos
    fase = 2*pi*n/NumSegmentos;
    
    %r_lab_vect = r_lab * [cos(fase), sin(fase), 0];
    v_lab_vect = v_lab * [-sin(fase), cos(fase), 0];
    %a_lab_vect = a_lab * [-cos(fase), -sin(fase), 0];
    
    % Velocidad del punto en esta fase.
    v_otro = Vel_Addition_Law(v_lab_vect, v_otro_SRI);
    
    % Ignoramos el 0, el 0 es tratado al final.
    if n > 0
               
        dv_otro = v_otro - v_otro_old;
        [eje_rot, ang] = ThomasPrecessionAng(Sigma_otro, dv_otro, v_otro);

        wT_otro = wT_otro + eje_rot * ang * (w_otro)/(2*pi);      
        
        %AxisRotMatrix = GetThomasRotMatrix(v_lab_vect, v_otro_SRI );
        %eje_rot = eje_rot * AxisRotMatrix;
        
        % El factor 2*pi es por no comparar un valor concreto contra una
        % longitud total que es 2*pi.
        %RotMatrix = RotationMatrixGeneral(eje_rot, ((ang*w_otro)/(2*pi)));
        
        % Con la versión original al menos no hacía locuras.
        %MatrizRotacionTotal = MatrizRotacionTotal * RotMatrix;
        
    end
    
    v_otro_old = v_otro;
    n = n + 1;
end

wT_otro = wT_otro
lT_otro2 = wT_otro * r_lab * Sigma_otro

%MatrizRotacionTotal = MatrizRotacionTotal

% Ya está bien que esta línea falle
COMPLETA DERROTA!

