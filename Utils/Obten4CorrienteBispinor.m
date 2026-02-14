function j = Obten4CorrienteBispinor(Phi)
    [jt, jx, jy, jz]  = ObtenCorrientesBispinor(Phi);
    
    j = [jt, jx, jy, jz];
end