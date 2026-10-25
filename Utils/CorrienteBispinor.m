function J_4 = CorrienteBispinor(bispinor)

[gt, gx, gy, gz] = MatricesGamma();

rho = bispinor * (gt * gt) * bispinor';

jx = bispinor * (gt * gx) * bispinor';

jy = bispinor * (gt * gy) * bispinor';

jz = bispinor * (gt * gz) * bispinor';

J_4 = [rho, jx, jy, jz];

end
