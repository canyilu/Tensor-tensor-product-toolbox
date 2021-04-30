function bX = bcirc(X)

% The block circulant matrix of a 3 way tensor
%
% X     -    n1*n2*n3 tensor
%
% bX    -    block circulant matrix of size (n1*n3)*(n2*n3)
%
% version 1.0 - 15/09/2017
%
% Written by Canyi Lu (canyilu@gmail.com)
% 
% We use the Corollary 4 in the paper
% Brian J Olson, Steven W Shaw, Chengzhi Shi, Christophe Pierre, Robert G Parker,
% Circulant matrices and their application to vibration analysis, Applied
% Mechanics Reviews, 2014
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
s = eye(n3,n3);
bX = zeros(n1*n3,n2*n3);
for i = 1 : n3
    S = gallery('circul',s(i,:)')';
    bX = bX + kron(S,X(:,:,i));
end