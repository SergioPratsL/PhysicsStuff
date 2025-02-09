function [r_ret, v_ret] = Obten_Posicion_Ret_O_Avz(r, v, a, Signo)
% La aceleración es la aceleración en el SRI propia

if Signo == 'a'
    factor = 1;
else
    factor = -1;
end


iteraciones = 20;
segmentos = 1000;

dist = norm(r);

n = 1;
dist_new = dist;

while n <= iteraciones    
    m = 1;
    t_seg = dist_new / segmentos;
    % Reiniciar las variables de velocidad y desplazamiento.
    v_t = v;
    despl = [0,0,0];
    
    % Actualizamos la velocidad y la posición a través de los segmentos
    % (integramos)
    while m <= segmentos
        despl = despl + v_t * t_seg * factor;
        
        Sigma = fGamma(v_t);
        p = v_t * Sigma;
        p = p + a * t_seg * factor;
        Sigma = fGammaFrom_p(p);
        v_t = p / Sigma;      
        
        m = m + 1;
    end
    
    % Signo menos porque es el origen del vector lo que se mueve!
    r_new = r - despl;
    dist_new = norm(r_new);
    n = n + 1;
end

r_ret = r_new;
v_ret = v_t;