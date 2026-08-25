
function R = MatrizRotacion(axis, angle_rad)
    % GET_ROTATION_MATRIX Computes the 3x3 spatial rotation matrix around an arbitrary axis.
    %
    % Inputs:
    %   axis      - 1x3 or 3x1 vector defining the spatial axis of rotation
    %   angle_rad - Rotation angle in radians (follows right-hand rule)
    %
    % Output:
    %   R         - 3x3 orthogonal rotation matrix (det(R) = 1)

    % Ensure input axis is a 3-element vector
    axis = axis(:);
    if length(axis) ~= 3
        error('Axis must be a 3-element vector.');
    end

    % Normalize the rotation axis vector
    axis_norm = norm(axis);
    if axis_norm == 0
        error('Rotation axis cannot be a zero vector.');
    end
    k = axis / axis_norm;

    % Convert angle to radians
    %theta = deg2rad(angle_deg);
    theta = angle_rad;

    % Construct the cross-product skew-symmetric matrix K
    K = [  0,   -k(3),  k(2);
         k(3),    0,   -k(1);
        -k(2),  k(1),    0  ];

    % Matrix form of Rodrigues' Rotation Formula: R = I + sin(theta)*K + (1 - cos(theta))*K^2
    R = eye(3) + sin(theta) * K + (1 - cos(theta)) * (K * K);
end