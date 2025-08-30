function A_prt = NuevoTerminoEM_Dirac_Prt(D_A, LCM)
% Las derivadas de los potenciales deben estar normalizadas para que c=1,
% de esta forma los campos E y B tienen las mismas unidades, de lo
% contrario a los términos magnéticos se les debería dividir por c*2
% mientras que a los electricos por c, lo cual lo complica todo..

    h = 6.626070 * 10^-34;
    h_bar = h / (2*pi);
    c = 299792458;
    m_e = 9.10938 * 10^-31;
    
    cof = h_bar / (m_e * c);

    V_prt = NuevoTerminoEM_index(1, D_A, LCM) * cof;
    Ax_prt = NuevoTerminoEM_index(2, D_A, LCM) * cof ;
    Ay_prt = NuevoTerminoEM_index(3, D_A, LCM) * cof;
    Az_prt = NuevoTerminoEM_index(4, D_A, LCM) * cof;

    A_prt = {V_prt, Ax_prt, Ay_prt, Az_prt};
end

function resultado = NuevoTerminoEM_index(index, D_A, LCM)
% index es la coordenada que estamos devolviendo y D_A la matriz de
% derivadas parciales del potencial

    gt = MatrizGamma(0);
    g5 = MatrizGamma(5);

    gi5  = gt * g5;

    m = 1;
    resultado = zeros(4);
    while m <= 4
        if (m == index)
            m = m + 1;
            continue;
        end

        n = 1;
        while n <= 4
            if (n == index || n == m)
                n = n + 1;
                continue;
            end

            o = 1;
            while o <= 4
                if (o == index || o == m || o == n)
                    o = o + 1;
                    continue
                end

                dpot = D_A(m,n);
                
                if( dpot ~= 0)
                    signo = LCM(index, o, m, n);
                    go = MatrizGamma((o-1));
                    resultado = resultado + signo * gi5 * go * dpot * factorMinkowski(m,n);
                end
                    
                o = o + 1;
            end

            n = n + 1;
        end

        m = m + 1;
    end
    
    % "Metrica de Minkowski"
    if index > 1
        %resultado = - resultado;
    end
    %resultado = (-1i) * resultado;
end

function factor = factorMinkowski(m,n)
    if m == 1
        factor1 = 1;
    else
        factor1 = -1;
    end

    if n == 1
        factor2 = 1;
    else
        factor2 = -1;
    end    
    
    factor = factor1 * factor2;
end
