function Xbdiag = bdiag(X)
%
% reformulate a 3 way tensor as a block diagonal matrix 
% X      - n1*n2*n3 tensor
% Xbdiag - (n1n3)*(n2n3) matrix
%
% version 1.0 - 18/06/2016
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

[n1,n2,n3] = size(X);
Xbdiag = zeros(n1*n3,n2*n3);

for i = 1 : n3
    Xbdiag((i-1)*n1+1:i*n1,(i-1)*n2+1:i*n2) = X(:,:,i);
end
