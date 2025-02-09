function Val = SumaBase( Vect, dif1, dif2)

Val_aux = Vect(1,1) * dif1(1) * dif2(1) + Vect(1,2) * dif2(1) * dif2(2) + Vect(1,3) * dif1(1) * dif2(3);

Val_aux = Val_aux + Vect(2,1) * dif1(2) * dif2(1) + Vect(2,2) * dif2(2) * dif2(2) + Vect(2,3) * dif1(2) * dif2(3);

Val_aux = Val_aux + Vect(3,1) * dif1(3) * dif2(1) + Vect(3,2) * dif2(3) * dif2(2) + Vect(3,3) * dif1(3) * dif2(3);

Val = Val_aux;
