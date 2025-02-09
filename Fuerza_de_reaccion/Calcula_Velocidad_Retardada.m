function v_ret = CalculaVelocidadRetardada(v, a, r, signo)
% Estos cálculos necesitan ser relativistas

if signo == 'r'
    v_ret = v - a * norm(rt_ret);
else
    v_avz = v + a * norm(rt_avz);
end


