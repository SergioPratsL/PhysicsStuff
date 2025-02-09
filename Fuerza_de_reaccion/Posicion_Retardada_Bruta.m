function rt = Posicion_Retardada_Bruta(r, v, a)

dist = norm(r);

n = 1;
dist_new = dist;
v_new = v;
while n < 100
    % El origen se atrasa dist_new, atrasar el origen igual a avanzar el
    % destino
    r_new = r + v_new * dist_new;
    dist_new = norm(r_new);
    v_new = v - (dist_new * a)/2;
    n = n + 1;
end

rt = r_new;