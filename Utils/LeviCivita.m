
% The Levi-Civita symbol, epsilon, the permutation tensor.
%--------------------------------------------------------------------------
% The rest of the file is the function itself.
% Bill Davidson, quellen@yahoo.com
% 13 Nov 2005
%
% Call it in this way using a vector of integers:
%   levi_civita_symbol([i j]), or
%   levi_civita_symbol([i j k]), or
%   levi_civita_symbol([i j k l]), etc ...
%
function [n]=LeviCivita(a) % a=[i j k ...]
if sum(sort(a)~=(1:length(a)))>0
    n=0;
else
    n=(-1)^mod(permutation_count(a),2);
end

% a permutation count
function n=permutation_count(a) % counts permutations
[b,n]=bubble_sort(a);
% a bubble sort with a permutation count

function [b,n]=bubble_sort(a) % n = a permutation count
b=a;
n=0;
m=length(a);
if m>1
    for i=1:m-1
        if b(i)>b(i+1)
            c=b(i);b(i)=b(i+1);b(i+1)=c;
            n=n+1;
        end
    end
    [c,k]=bubble_sort(b(1:m-1));
    b=[c b(m)];
    n=n+k;
end
