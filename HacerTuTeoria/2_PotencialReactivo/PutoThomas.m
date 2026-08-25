function [vBC_rot] = PutoThomas(vBA, vBC)
% Borradme por favor :P.

    vCA = Vel_Addition_Law( vBA, vBC );
    vAC = Vel_Addition_Law( vBC, vBA );
    eje_rot = cross(vBC, vBA);
    eje_rot = eje_rot / norm(eje_rot);    
    val = -dot(vAC, vCA) / (vCA(1)^2+vCA(2)^2+vCA(3)^2);
    % Sin esto, errores de precisión podían devolver un valor imaginario
    if( abs(val) > 1)
        val = sign(val);
    end
    angulo = acos(val);
    
    % Sólo funciona en ejes X, Y, quizá no sirva para nada...
    vBC_rot(1) = vBC(1) * cos(angulo) + vBC(2) * sin(angulo);
    vBC_rot(2) = vBC(2) * cos(angulo) - vBC(1) * sin(angulo);
    vBC_rot(3) = vBC(3);
    
    