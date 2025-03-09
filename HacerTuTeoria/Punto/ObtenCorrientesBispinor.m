function [jt, jx, jy, jz] = ObtenCorrientesBispinor(Phi_ori)
    
    Phi = Phi_ori.';

    [alfa_t, alfa_x, alfa_y, alfa_z] = MatricesAlfa();
    
    jt = Phi' * alfa_t * Phi;
    jx = Phi' * alfa_x * Phi;
    jy = Phi' * alfa_y * Phi;
    jz = Phi' * alfa_z * Phi;    

end