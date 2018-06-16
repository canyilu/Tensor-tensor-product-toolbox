function tsn = tsn(X)

% Tensor spectral norm of a 3 way tensor
%
% X     - n1*n2*n3 tensor
% tsn    - tensor spectral norm
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

n3 = size(X,3);
X = fft(X,[],3);

% i=1
tsn = norm(X(:,:,1),2);

% i=2,...,halfn3
halfn3 = round(n3/2);
for i = 2 : halfn3
    tsn = max(tsn,norm(X(:,:,i),2));
end

% if n3 is even
if mod(n3,2) == 0
    i = halfn3+1;
    tsn = max(tsn,norm(X(:,:,i),2));
end