function invX = tinv(X)

% tinv(X) is the inverse of the tensor X of size n*n*n3.
%   A warning message is printed if X is badly scaled or
%   nearly singular.
%
% version 1.0 - 14/06/2018
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
% Norm, arXiv preprint arXiv:1804.03728, 2018
%

[n1,n2,n3] = size(X);
if n1 ~= n2
    error('Error using tinv. Tensor must be square.');
end

X = fft(X,[],3);
invX = zeros(n1,n2,n3);
I = eye(n1);
% first frontal slice
invX(:,:,1) = X(:,:,1)\I;
% i=2,...,halfn3
halfn3 = round(n3/2);
for i = 2 : halfn3
    invX(:,:,i) = X(:,:,i)\I;
    invX(:,:,n3+2-i) = conj(invX(:,:,i));
end
% if n3 is even
if mod(n3,2) == 0
    i = halfn3+1;
    invX(:,:,i) = X(:,:,i)\I;
end
invX = ifft(invX,[],3);

