
clear

syms rx;
syms ry;
syms rz;

syms vx;
syms vy;
syms vz;


  
V = 1 / sqrt(rx^2+ry^2+rz^2);
 
Ax1 = 1/sqrt(rx^2+ry^2+rz^2) * ( vx*(rx^2+ry^2+rz^2) - rx*(rx*vx+ry*vy+rz*vz) ) / ((rx^2+ry^2+rz^2)*(vx^2+vy^2+vz^2)-(rx*vx+ry*vy+rz*vz)^2);

Ay1 = 1/sqrt(rx^2+ry^2+rz^2) * ( vy*(rx^2+ry^2+rz^2) - ry*(rx*vx+ry*vy+rz*vz) ) / ((rx^2+ry^2+rz^2)*(vx^2+vy^2+vz^2)-(rx*vx+ry*vy+rz*vz)^2);

Az1 = 1/sqrt(rx^2+ry^2+rz^2) * ( vz*(rx^2+ry^2+rz^2) - ry*(rx*vx+ry*vy+rz*vz) ) / ((rx^2+ry^2+rz^2)*(vx^2+vy^2+vz^2)-(rx*vx+ry*vy+rz*vz)^2);



Ax2 = 1/sqrt(vx^2+vy^2+vz^2) * ( rx*(vx^2+vy^2+vz^2) - vx*(rx*vx+ry*vy+rz*vz)) / ( (vy*rz-vz*ry)^2 + (vz*rx-vx*rz)^2 + (vx*ry-vy*rx)^2 );

Ay2 = 1/sqrt(vx^2+vy^2+vz^2) * ( ry*(vx^2+vy^2+vz^2) - vy*(rx*vx+ry*vy+rz*vz)) / ( (vy*rz-vz*ry)^2 + (vz*rx-vx*rz)^2 + (vx*ry-vy*rx)^2 );

Az2 = 1/sqrt(vx^2+vy^2+vz^2) * ( rz*(vx^2+vy^2+vz^2) - vz*(rx*vx+ry*vy+rz*vz)) / ( (vy*rz-vz*ry)^2 + (vz*rx-vx*rz)^2 + (vx*ry-vy*rx)^2 );



Ax = Ax1 + Ax2;

Ay = Ay1 + Ay2;

Az = Az1 + Az2;

