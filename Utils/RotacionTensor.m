function T_rot = RotacionTensor(T, dir_rot, angulo)

    mat_rot_4_4 = MatrizRotacion_4_4(dir_rot, angulo);
    
    T_rot = mat_rot_4_4 * T * mat_rot_4_4';
end