function [axis, angle] = ExtraeAnguloDeMatriz(R_ini)
    % ROTM2AXISANGLE Extracts the rotation axis and angle from a 3x3 rotation matrix.
    %   [AXIS, ANGLE] = ROTM2AXISANGLE(R) takes a 3x3 rotation matrix R (SO(3))
    %   and returns:
    %     AXIS  - 3x1 normalized column vector representing the axis of rotation.
    %     ANGLE - Scalar rotation angle in radians, in the range [0, pi].


    if isequal(size(R_ini), [4, 4])
        R = R_ini(2:4, 2:4);
    else
        R = R_ini;
    end
    
    % Input validation
    if ~isequal(size(R), [3, 3])
        error('Input must be a 3x3 matrix.');
    end

    % 1. Compute the rotation angle from the trace of R
    %    trace(R) = 1 + 2*cos(theta)
    tr = trace(R);
    cos_theta = (tr - 1) / 2;
    
    % Clamp cos_theta to [-1, 1] to guard against floating-point inaccuracies
    cos_theta = max(-1, min(1, cos_theta));
    angle = acos(cos_theta);

    % 2. Extract the rotation axis based on the angle
    if angle < 1e-6
        % Case A: Angle is approximately 0 (Identity matrix, no rotation)
        % Axis can be any unit vector; default to [1; 0; 0]
        axis = [1; 0; 0];
        angle = 0;
        
    elseif abs(angle - pi) < 1e-6
        % Case B: Angle is approximately pi (180 degrees)
        % R = 2*u*u' - I, so u*u' = (R + I) / 2
        % Find the column of (R + I) with the largest norm for numerical stability
        M = R + eye(3);
        [~, max_col] = max(sum(M.^2, 1));
        axis = M(:, max_col);
        axis = axis / norm(axis);
        angle = pi;
        
    else
        % Case C: Standard case (0 < angle < pi)
        % The skew-symmetric part gives: [rx; ry; rz] = 1 / (2*sin(theta)) * [R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)]
        v = [R(3,2) - R(2,3);
             R(1,3) - R(3,1);
             R(2,1) - R(1,2)];
        
        axis = v / norm(v);
    end
end