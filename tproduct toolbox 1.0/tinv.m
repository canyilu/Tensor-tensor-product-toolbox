function X = tinv(X)

% tinv(X) is the inverse of the tensor X of size n*n*n3.
%   A warning message is printed if X is badly scaled or
%   nearly singular.
%
% version 1.0 - 14/06/2018
% version 1.1 - 28/04/2021
%
% Written by Canyi Lu (canyilu@gmail.com)
%
% References: 
% Canyi Lu, Tensor-Tensor Product Toolbox. Carnegie Mellon University. 
% June, 2018. https://github.com/canyilu/tproduct.
%
% Canyi Lu, Jiashi Feng, Yudong Chen, Wei Liu, Zhouchen Lin and Shuicheng
% Yan, Tensor Robust Principal Component Analysis with A New Tensor Nuclear
% Norm, TPAMI, 2019
%

[n1,n2,n3] = size(X);
if n1 ~= n2
    error('Error using tinv. Tensor must be square.');
end

X = fft(X,[],3);
I = eye(n1);
halfn3 = ceil((n3+1)/2);
for i = 1 : halfn3
    X(:,:,i) = X(:,:,i)\I;
end
for i = halfn3+1 : n3
    X(:,:,i) = conj(X(:,:,n3+2-i));
end
X = ifft(X,[],3);
    
