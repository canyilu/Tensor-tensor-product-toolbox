function C = tprod(A,B)

% Tensor-tensor product of two 3 way tensors: C = A*B
% A - n1*n2*n3 tensor
% B - n2*l*n3  tensor
% C - n1*l*n3  tensor
%
% version 2.0 - 09/10/2017
% version 2.1 - 28/04/2021
%
% Written by Canyi Lu (canyilu@gmail.com)
%
%
% References: 
% Canyi Lu, Tensor-Tensor Product Toolbox. Carnegie Mellon University. 
% June, 2018. https://github.com/canyilu/tproduct.
%
% Canyi Lu, Jiashi Feng, Yudong Chen, Wei Liu, Zhouchen Lin and Shuicheng
% Yan, Tensor Robust Principal Component Analysis with A New Tensor Nuclear
% Norm, TPAMI, 2019
%

[n1,n2,n3] = size(A);
[m1,m2,m3] = size(B);

if n2 ~= m1 || n3 ~= m3 
    error('Inner tensor dimensions must agree.');
end

A = fft(A,[],3);
B = fft(B,[],3);
C = zeros(n1,m2,n3);

halfn3 = ceil((n3+1)/2);
for i = 1 : halfn3
    C(:,:,i) = A(:,:,i)*B(:,:,i);
end
for i = halfn3+1 : n3
    C(:,:,i) = conj(C(:,:,n3+2-i));
end
C = ifft(C,[],3);
