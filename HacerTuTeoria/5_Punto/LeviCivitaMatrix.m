function L_C = LeviCivitaMatrix(dims) 
% Matlab parece no saber tratar muy bien las matrices de rango más alto,
% saco esto de Utils porque más bien es IN-util.

if (dims > 4 || dims < 2)
    disp('El número de dimensiones debe estar entre 2 y 4');
end

if dims == 2
    L_C = leviCivitaMatrixDim2();
end

if dims == 3
    L_C = leviCivitaMatrixDim3();
end

if dims == 4
    L_C = leviCivitaMatrixDim4();
end


function L_C = leviCivitaMatrixDim2()

L_C = [0, 1; -1, 0];


function L_C = leviCivitaMatrixDim3()

L_C = zeros(3,3,3);

for i=1:3
    for j=1:3
        if j == i
            continue
        end
        for k=1:3
            if k == j || k == i
                continue
            end
           
            L_C(i,j,k) = LeviCivita([i,j,k]);
        end
    end
end


function L_C = leviCivitaMatrixDim4()

L_C = zeros(4,4,4,4);

for i=1:4
    for j=1:4
        if j == i
            continue
        end
        for k=1:4
            if k == j || k == i
                continue
            end
            
            for m=1:4            
                if m == k || m == j || m == i
                    continue
                end                
                L_C(i,j,k,m) = LeviCivita([i,j,k,m]);
            end
        end
    end
end
