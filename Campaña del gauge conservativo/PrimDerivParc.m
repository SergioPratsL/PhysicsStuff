function [d1_V, d1_Ax, d1_Ay, d1_Az] = PrimDerivParc(Rx, Ry, v)

DPot_x = DerivadasPotencialYinn( Rx, Ry, v, 'x');
DPot_y = DerivadasPotencialYinn( Rx, Ry, v, 'y');
DPot_z = DerivadasPotencialYinn( Rx, Ry, v, 'z');

d1_V = [DPot_x(1), DPot_y(1), DPot_z(1)];

d1_Ax = [DPot_x(2), DPot_y(2), DPot_z(2)];

d1_Ay = [DPot_x(3), DPot_y(3), DPot_z(3)];

d1_Az = [DPot_x(4), DPot_y(4), DPot_z(4)];



