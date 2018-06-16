function tnn = tnn(X)

% Tensor nuclear norm of a 3 way tensor
%
% X     - n1*n2*n3 tensor
% tnn   - tensor nuclear norm
%
% version 2.0 - 09/10/2017
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

n3 = size(X,3);
X = fft(X,[],3);
tnn = 0;

% i=1
s = svd(X(:,:,1),'econ');
tnn = tnn+sum(s);

% i=2,...,halfn3
halfn3 = round(n3/2);
for i = 2 : halfn3
    s = svd(X(:,:,i),'econ');
    tnn = tnn+sum(s)*2;
end

% if n3 is even
if mod(n3,2) == 0
    i = halfn3+1;
    s = svd(X(:,:,i),'econ');
    tnn = tnn+sum(s);
end
tnn = tnn/n3;
