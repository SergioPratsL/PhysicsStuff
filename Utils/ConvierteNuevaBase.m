function vector_nb = ConvierteNuevaBase( vector, e1, e2, e3 )
% Las bases e1, e2 y e3 deben estar normalizadas!!

vector_nb = [dot(vector,e1), dot(vector,e2), dot(vector,e3)];

end