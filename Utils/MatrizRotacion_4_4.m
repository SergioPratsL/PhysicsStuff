
function R4 = MatrizRotacion_4_4(axis, angle_rad)
    % GET_LORENTZ_ROTATION_MATRIX Returns a 4x4 spacetime transformation matrix 
    % for a pure spatial rotation around an arbitrary axis.
    %
    % Inputs:
    %   axis      - 1x3 or 3x1 vector defining the spatial axis of rotation
    %   angle_rad - Rotation angle in radians (follows right-hand rule)
    %
    % Output:
    %   R4        - 4x4 Lorentz rotation matrix acting on 4-vectors [c*t; x; y; z]

    % Get the 3x3 spatial rotation matrix
    R3 = MatrizRotacion(axis, angle_rad);

    % Construct the 4x4 spacetime rotation matrix
    R4 = eye(4);
    R4(2:4, 2:4) = R3;
end