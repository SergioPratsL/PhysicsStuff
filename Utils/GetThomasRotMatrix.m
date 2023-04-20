function [Matrix3, angulo] = GetThomasRotMatrix(vBA, vBC)

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
    Matrix3 = RotationMatrixGeneral(eje_rot, angulo);
    