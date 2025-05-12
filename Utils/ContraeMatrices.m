
function val = ContraeMatrices(T1, T2)

    dim_T1 = size(T1);
    dim_T2 = size(T1);

    if(dim_T1(1) ~= dim_T1(2) || dim_T1(1) ~= dim_T2(1) || dim_T1(1) ~= dim_T2(2))
        error('Las matrices deben ser cuadradas y de la misma dimensión')
    end

    dimensiones = dim_T1(1);

    val = 0;

    m = 1;
    while m <= dimensiones
        n = 1;
        while n <= dimensiones
            val = val + T1(m,n) * T2(m,n);
            n = n + 1;
        end

        m = m + 1;
    end
