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
X = fft(X,[],3);

if n1>n2 && exist('opt', 'var') && strcmp(opt,'econ') == 1    
    Q = zeros(n1,n2,n3);
    R = zeros(n2,n2,n3);    
    % first frontal slice
    [Q(:,:,1),R(:,:,1)] = qr(X(:,:,1),0);
    % i=2,...,halfn3
    halfn3 = round(n3/2);
    for i = 2 : halfn3
        [Q(:,:,i),R(:,:,i)] = qr(X(:,:,i),0);
        Q(:,:,n3+2-i) = conj(Q(:,:,i));
        R(:,:,n3+2-i) = conj(R(:,:,i));
    end
    % if n3 is even
    if mod(n3,2) == 0
        i = halfn3+1;
        [Q(:,:,i),R(:,:,i)] = qr(X(:,:,i),0);
    end
    
else    
    Q = zeros(n1,n1,n3);
    R = zeros(n1,n2,n3);
    % first frontal slice
    [Q(:,:,1),R(:,:,1)] = qr(X(:,:,1));
    % i=2,...,halfn3
    halfn3 = round(n3/2);
    for i = 2 : halfn3
        [Q(:,:,i),R(:,:,i)] = qr(X(:,:,i));
        Q(:,:,n3+2-i) = conj(Q(:,:,i));
        R(:,:,n3+2-i) = conj(R(:,:,i));
    end
    % if n3 is even
    if mod(n3,2) == 0
        i = halfn3+1;
        [Q(:,:,i),R(:,:,i)] = qr(X(:,:,i));
    end
end

Q = ifft(Q,[],3);
R = ifft(R,[],3);
