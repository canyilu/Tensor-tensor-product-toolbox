function Xbdiag = bdiag(X)
%
% reformulate a 3 way tensor as a block diagonal matrix 
% X      - n1*n2*n3 tensor
% Xbdiag - (n1n3)*(n2n3) matrix
%
% version 1.0 - 18/06/2016
% version 1.1 - 07/11/2018
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


sizen = size(X);
Xbdiag = zeros(sizen(1)*sizen(3),sizen(2)*sizen(3));

for i = 1 : sizen(3)
    Xbdiag((i-1)*sizen(1)+1:i*sizen(1),(i-1)*sizen(2)+1:i*sizen(2)) = X(:,:,i);
end
