function f_new = AproximacionTaylor5( f, dtf, dt2f, dt3f, dt4f, dt )
% Devuelve la aproximación de Taylor usando hasta la cuarta derivada

dt2 = dt^2;
dt3 = dt2 * dt;
dt4 = dt3 * dt;

f_new = f + dtf * dt + 0.5 * dt2f * dt2 + (1/6) * dt3f * dt3 + (1/24) * dt4f * dt4;

end

