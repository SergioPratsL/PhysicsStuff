function v = DeterminaSric(E, B)
% DeterminaSric calculates the Lorentz boost velocity required 
% to reach a reference frame where the Poynting vector is zero.
%
% Inputs:
%   E - Electric field vector [Ex, Ey, Ez] in happy units
%   B - Magnetic field vector [Bx, By, Bz] in happy units
%
% Output:
%   v - Velocity vector [vx, vy, vz] normalized to speed of light

    % Ensure inputs are 3-element vectors
    if numel(E) ~= 3 || numel(B) ~= 3
        error('Inputs E and B must be 3-element vectors.');
    end
    
    % Store original shape to return output in the same format
    input_size = size(E);
    
    % Convert to column vectors for uniform calculation
    E = E(:);
    B = B(:);
        
    % Calculate the cross product and its magnitude
    ExB = cross(E, B);
    S_mag = norm(ExB);
    
    % If the Poynting vector is already zero, v = 0
    if S_mag == 0
        v = zeros(input_size);
        return;
    end
    
    % Calculate the magnitudes squared of the fields
    E_mag2 = sum(E.^2);
    B_mag2 = sum(B.^2);
    
    % Calculate the dimensionless parameter alpha
    alpha = (S_mag) / (E_mag2 + B_mag2);
    
    % Handle the plane-wave limit (alpha ~ 0.5)
    if alpha >= (0.5 - 1e-8)
        warning('Electromagnetic plane wave detected (E = cB, E perp B). Required velocity is c.');
        beta = 1;
    else
        % Solve for beta (v/c)
        beta = (1 - sqrt(1 - 4 * alpha^2)) / (2 * alpha);
    end
    
    % Calculate final velocity vector
    v_mag = beta;
    v_dir = ExB / S_mag;
    v_vec = v_mag * v_dir;
    
    % Reshape output to match the dimension format of the input E
    v = reshape(v_vec, input_size);
end
