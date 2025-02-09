function [ E, B ] = ObtenCampoEMTensor( T )
%OBTENCAMPOEMTENSOR Summary of this function goes here
%   Detailed explanation goes here

E = [T(2,1), T(3,1), T(4,1)];

B = [T(3,4), T(4,2), T(2,3)];

end

