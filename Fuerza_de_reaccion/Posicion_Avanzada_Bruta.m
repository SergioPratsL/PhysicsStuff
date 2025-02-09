function rt = Posicion_Avanzada_Bruta(r, v, a)

dist = norm(r);

n = 1;
dist_new = dist;
v_new = v;
while n < 100
    r_new = r - v_new * dist_new;
    dist_new = norm(r_new);
    v_new = v + (dist_new * a)/2;
    n = n + 1;
end

rt = r_new;