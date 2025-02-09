function [alfa_t, alfa_x, alfa_y, alfa_z] = MatricesAlfa()
% Devuelve las 4 matrices gamma, ya que generalmente las necesito todas.

[gt, gx, gy, gz] = MatricesGamma();

alfa_t = gt * gt;
alfa_x = gt * gx;
alfa_y = gt * gy;
alfa_z = gt * gz;