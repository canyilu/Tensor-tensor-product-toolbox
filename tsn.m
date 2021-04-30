function tsn = tsn(X)

% Tensor spectral norm of a 3 way tensor
%
% X     - n1*n2*n3 tensor
% tsn   - tensor spectral norm
%
% version 1.0 - 14/06/2018
% version 1.1 - 28/04/2021 a more efficient version
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

X = fft(X,[],3);
tsn = 0;
for i = 1 : ceil((size(X,3)+1)/2)
    tsn = max(tsn,norm(X(:,:,i),2));
end

