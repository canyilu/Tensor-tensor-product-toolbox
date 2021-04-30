function [X, tnn] = prox_tnn(Y,rho,transform)

% Tensor singular value thresholding (TSVT) which solves the following
% problem (the proximal operator of the tensor nuclear norm under linear transform)
%
% min_X rho*||X||_*+0.5*||X-Y||_F^2
%
% An error message is printed if the property L'*L=L*L'=l*I, for some l>0,
% of the linear transform does not hold.
%
%
% Input:
%       Y       -   n1*n2*n3 tensor
%       rho     -   a constant > 0
%   transform   -   a structure which defines the linear transform
%       transform.L: the linear transform of two types:
%                  - type I: function handle, i.e., @fft, @dct
%                  - type II: invertible matrix of size n3*n3
%
%       transform.inverseL: the inverse linear transform of transform.L
%                         - type I: function handle, i.e., @ifft, @idct
%                         - type II: inverse matrix of transform.L
%
%       transform.l: a constant which indicates whether the following property holds for the linear transform or not:
%                    L'*L=L*L'=l*I, for some l>0.                           
%                  - transform.l > 0: indicates that the above property holds. Then we set transform.l = l.
%                  - transform.l < 0: indicates that the above property does not hold. Then we can set transform.l = c, for any constant c < 0.
%       If not specified, fft is the default transform, i.e.,
%       transform.L = @fft, transform.l = n3, transform.inverseL = @ifft. 
%
%
% Output:
%        X      -   the solution of tensor singular value thresholding
%        tnn    -   tensor nuclear norm of X
%
%
%
% See also tnn, tsvd, lineartransform, inverselineartransform
%
%
% References:
% Canyi Lu, Tensor-Tensor Product Toolbox. Carnegie Mellon University. 
% June, 2018. https://github.com/canyilu/tproduct.
%
% Canyi Lu, Xi Peng, Yunchao Wei, Low-Rank Tensor Completion With a New Tensor 
% Nuclear Norm Induced by Invertible Linear Transforms. IEEE International 
% Conference on Computer Vision and Pattern Recognition (CVPR), 2019
%
% Canyi Lu, Pan Zhou. Exact Recovery of Tensor Robust Principal Component 
% Analysis under Linear Transforms. arXiv preprint arXiv:1907.08288. 2019
%
%
% version 1.0 - 01/02/2019
% version 1.1 - 29/04/2021
%
% Written by Canyi Lu (canyilu@gmail.com)
%

[n1,n2,n3] = size(Y);
if nargin == 3
    if transform.l < 0
        error("property L'*L=L*L'=l*I does not holds for some l>0.");
    end
else    
    % fft is the default transform
    transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
end

X = zeros(n1,n2,n3);
tnn = 0;
if isequal(transform.L,@fft)
    % efficient computing for fft transform
    Y = fft(Y,[],3);    
    % first frontal slice
    [U,S,V] = svd(Y(:,:,1),'econ');
    S = diag(S);
    r = length(find(S>rho));
    if r >= 1
        S = S(1:r)-rho;
        X(:,:,1) = U(:,1:r)*diag(S)*V(:,1:r)';
    end
    tnn = tnn+sum(S);
    % i=2,...,halfn3
    halfn3 = round(n3/2);
    for i = 2 : halfn3
        [U,S,V] = svd(Y(:,:,i),'econ');
        S = diag(S);
        r = length(find(S>rho));
        if r >= 1
            S = S(1:r)-rho;
            X(:,:,i) = U(:,1:r)*diag(S)*V(:,1:r)';
        end
        X(:,:,n3+2-i) = conj(X(:,:,i));
        tnn = tnn+sum(S)*2;
    end
    % if n3 is even
    if mod(n3,2) == 0
        i = halfn3+1;
        [U,S,V] = svd(Y(:,:,i),'econ');
        S = diag(S);
        r = length(find(S>rho));
        if r >= 1
            S = S(1:r)-rho;
            X(:,:,i) = U(:,1:r)*diag(S)*V(:,1:r)';
        end
        tnn = tnn+sum(S);
    end
    tnn = tnn/n3;
    X = ifft(X,[],3);
else
    % other transform
    Y = lineartransform(Y,transform);
    for i = 1 : n3
        [U,S,V] = svd(Y(:,:,i),'econ');
        S = diag(S);
        r = length(find(S>rho));
        if r >= 1
            S = S(1:r)-rho;
            X(:,:,i) = U(:,1:r)*diag(S)*V(:,1:r)';
            tnn = tnn+sum(S);
        end
    end
    tnn = tnn/transform.l;
    X = inverselineartransform(X,transform);
end
