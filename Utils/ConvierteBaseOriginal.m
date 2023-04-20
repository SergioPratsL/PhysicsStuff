function vector_ori = ConvierteBaseOriginal( vector, e1, e2, e3 )
% Las bases e1, e2 y e3 deben estar normalizadas!!

vector_ori_1 = [vector(1)*e1(1) + vector(2)*e2(1) + vector(3)*e3(1)];
vector_ori_2 = [vector(1)*e1(2) + vector(2)*e2(2) + vector(3)*e3(2)];
vector_ori_3 = [vector(1)*e1(3) + vector(2)*e2(3) + vector(3)*e3(3)];

vector_ori = [vector_ori_1, vector_ori_2, vector_ori_3];

end
