function [Q,R] = tqr(X,opt)

% Tensor orthogonal-triangular decomposition.
%   [Q,R] = tqr(X), where X is n1*n2*n3, produces an n1*n2*n3 upper triangular
%   tensor R and an n1*n1*n3 orthogonal tensor Q so that X = Q*R.
%
%   [Q,R] = tqr(X,'econ') produces the "economy size" decomposition.
%   If n1>n2, only the first n2 lateral slices of Q and the first n2 horizontal
%   slices of R are computed. If n1<=n2, this is the same as [Q,R] = tqr(X).
%
% version 1.0 - 14/06/2018
% version 1.1 - 28/04/2021
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
X = fft(X,[],3);

if n1>n2 && exist('opt', 'var') && strcmp(opt,'econ') == 1    
    Q = zeros(n1,n2,n3);
    R = zeros(n2,n2,n3);
    halfn3 = ceil((n3+1)/2);
    for i = 1 : halfn3
        [Q(:,:,i),R(:,:,i)] = qr(X(:,:,i),0);
    end
    for i = halfn3+1 : n3
        Q(:,:,i) = conj(Q(:,:,n3+2-i));
        R(:,:,i) = conj(R(:,:,n3+2-i));
    end
else    
    Q = zeros(n1,n1,n3);
    R = zeros(n1,n2,n3);
    halfn3 = ceil((n3+1)/2);
    for i = 1 : halfn3
        [Q(:,:,i),R(:,:,i)] = qr(X(:,:,i));
    end
    for i = halfn3+1 : n3
        Q(:,:,i) = conj(Q(:,:,n3+2-i));
        R(:,:,i) = conj(R(:,:,n3+2-i));
    end
end

Q = ifft(Q,[],3);
R = ifft(R,[],3);
