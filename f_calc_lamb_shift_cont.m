% Continuacion de f_calc_lamb_shit para poder debugar mejor
% No ejecutar individualmente!

% Hacemos las FFTs...

if sacar_fft == 1
    f_onda_p = fft( f_onda );

    clear f_onda;

% Conjugar la funcion de onda!!
    f_onda_p = conj(f_onda_p);

    f_onda_pot_p = fft( f_onda_pot ); 

    clear f_onda_pot;
    
% Matrices f_onda y f_onda_pot terminadas
    fecha_hora = clock;

    hora = round(10 * fecha_hora(4) / 10);
    minuto = round(10 * fecha_hora(5) / 10);
    segundo = round(10 *fecha_hora(6) / 10);

    fprintf('FFTs terminadas: %d:%d:%d \r\n', hora, minuto, segundo);
    
end


f_unit = c_planck_reduc * 1 / (r_unit * num_iter_cartesiana);
f_unit_cuadr = f_unit^2;

% f_unit_cuadr_c_cuadr = f_unit_cuadr * c_cuadr;    % Esto es incorrecto. 


% Variables de chequeo
Onda_p_cuadr_acum = 0;
Onda_V_p_cuadr_acum = 0;

max_valor = 0;
Valor_acumulado = 0;
cont_x = 1;
while cont_x <= num_iter_cartesiana
    cont_y = 1;
    
    while cont_y <= num_iter_cartesiana
        cont_z = 1;    
        
        while cont_z <= num_iter_cartesiana;
      
% Ojo! En la DFT las frecuencias por encima de la mitad son en verdad
% frecuencias negativas!
            if cont_x <= num_iter_cartesiana_mitad
                cont_x_aux = cont_x - 1;                
            else
                cont_x_aux = cont_x - num_iter_cartesiana - 1;   
            end

            if cont_y <= num_iter_cartesiana_mitad
                cont_y_aux = cont_y - 1;                
            else
                cont_y_aux = cont_y - num_iter_cartesiana - 1;   
            end            
            
            if cont_z <= num_iter_cartesiana_mitad
                cont_z_aux = cont_z - 1;                
            else
                cont_z_aux = cont_z - num_iter_cartesiana - 1;   
            end            
            
            
            momentum_cuadr = f_unit_cuadr * ( (cont_x_aux)^2 + (cont_y_aux)^2 + (cont_z_aux)^2);      
            Propagador = 1 / ( term_ener_cuadr + momentum_cuadr);               
    
            Onda_p_cuadr_acum = Onda_p_cuadr_acum + f_onda_p( cont_x, cont_y, cont_z ) * conj(f_onda_p( cont_x, cont_y, cont_z ));
            Onda_V_p_cuadr_acum = Onda_V_p_cuadr_acum + f_onda_pot_p( cont_x, cont_y, cont_z ) * conj(f_onda_pot_p( cont_x, cont_y, cont_z ));

                
            Valor = f_onda_p( cont_x, cont_y, cont_z ) * f_onda_pot_p( cont_x, cont_y, cont_z ) * Propagador;

            if Valor > max_valor;
                max_valor = Valor;
                if max_valor > 7.719e+054
                    max_valor = max_valor;
                end
            end
            
            Valor_acumulado = Valor_acumulado + Valor;
        
            cont_z = cont_z + 1;
        end
        
        cont_y = cont_y + 1;        
    end
    
    cont_x = cont_x + 1;
end


% Todo el proceso terminado
fecha_hora = clock;

hora = round(10 * fecha_hora(4) / 10);
minuto = round(10 * fecha_hora(5) / 10);
segundo = round(10 *fecha_hora(6) / 10);

fprintf('Fin del proceso: %d:%d:%d \r\n', hora, minuto, segundo);


% Aplicamos el factor de volumen de cada punto sobre el que hemos integrado
Valor_acumulado = Valor_acumulado * vol_element * e;            % No nos olvidemos de la e!!!

% Para que la normalizacion de igual no debe elevarse al cuadrado!
% Valor_acumulado = Valor_acumulado / num_iter_cartesiana^2;
% Valor_acumulado = Valor_acumulado / num_iter_cartesiana;


Lamb_Shift_Pratiano = Valor_acumulado       % :)

Lamb_Shift_Pratiano_frec = Lamb_Shift_Pratiano / h


Onda_p_cuadr_acum
Onda_cuadr_acum

Ratio = Onda_p_cuadr_acum / Onda_cuadr_acum

Onda_V_p_cuadr_acum
Onda_V_cuadr_acum

Norm_Check = Onda_p_cuadr_acum * vol_element

Norm_Check2 = Onda_cuadr_acum * vol_element


max_valor

sacar_fft = 0;