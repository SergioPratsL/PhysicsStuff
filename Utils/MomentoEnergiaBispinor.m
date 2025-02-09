function mom4 = MomentoEnergiaBispinor(bispinor)

[gt, gx, gy, gz] = MatricesGamma();

E = bispinor * (gt * gt) * bispinor';

px = bispinor * (gt * gx) * bispinor';

py = bispinor * (gt * gy) * bispinor';

pz = bispinor * (gt * gz) * bispinor';

mom4 = [E, px, py, pz];

end
