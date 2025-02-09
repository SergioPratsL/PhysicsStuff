% Prueabas locas, las hice lanzando los comandos a mano y todas salieron mal...

b = MatrizGamma(0);
Gx = MatrizGamma(1);
Gy = MatrizGamma(2);
Gz = MatrizGamma(3);
G5 = MatrizGamma(5);

M0 = [0,0;0,0];

px = MatrizPauli(1);

bpx = [px, M0; M0, px];

bpx_inv = [M0, px; px, M0];


x = MatrixConmut(b, bpx_inv);
% if (x == bpx || x == -bpx)
%     disp('1')
% end

x = MatrixConmut(b, bpx);
% if (x == bpx || x == -bpx)
%     disp('2')
% end

x = MatrixConmut(b, Gx);
% if (x == bpx || x == -bpx)
%     disp('3')
% end

x = MatrixConmut(b, G5);
% if (x == bpx || x == -bpx)
%     disp('4')
% end

x = MatrixConmut(b, Gx*G5);
% if (x == bpx || x == -bpx)
%     disp('5')
% end

x = MatrixConmut(b, Gx*bpx);
% if (x == bpx || x == -bpx)
%     disp('6')
% end

x = MatrixConmut(b, Gx*bpx_inv);
% if (x == bpx || x == -bpx)
%     disp('7')
% end
