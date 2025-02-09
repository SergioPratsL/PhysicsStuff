function exp_M = ExponencialMatriz(M, num_iter)

if (num_iter < 10)
    num_iter = 10;
end

dims = size(M);

if (dims(1) ~= dims(2))
    error('ExponencialMatriz: matrix dimensions must agree')
end

n = 0;
fact = 1;

M_pot = eye(dims(1));
exp_M = zeros(dims(1));

while n <= num_iter
    
    if n > 1
        fact = fact * n;
    end
    
    exp_M = exp_M + M_pot / fact;
    
    M_pot = M_pot * M;
    
    n = n + 1;
end
    

