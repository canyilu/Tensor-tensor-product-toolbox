function [U,S,V] = tsvd(A,transform,opt)

% [U,S,V] = tsvd(A,transform,opt) computes the tensor SVD under linear transform, i.e., A=U*S*V^*, where S
% is a f-diagonal tensor, U and V are orthogonal under linear transform.
%
%
% Input:
%       A       -   n1*n2*n3 tensor
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
%       opt     -   options for different outputs of U, S and V:
%                   'full': (default) produces full tensor SVD, i.e., A = U*S*V^*, where
%                       U - n1*n1*n3
%                       S - n1*n2*n3
%                       V - n2*n2*n3
%                   'econ': produces the "economy size" decomposition. 
%                       Let m = min(n1,n2). Then, A = U*S*V^*, where
%                       U - n1*m*n3
%                       S - m*m*n3
%                       V - n2*m*n3
%                   'skinny': produces the skinny tensor SVD.
%                       Let r be the tensor tubal rank of A. Then, A = U*S*V^*, where
%                       U - n1*r*n3
%                       S - r*r*n3
%                       V - n2*r*n3
%
%
% Output: U, S, V
%
%
%
% See also lineartransform, inverselineartransform
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

if nargin < 3
    opt = 'full';
end
[n1,n2,n3] = size(A);
if nargin < 2
    % fft is the default transform
    transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
end

if isequal(transform.L,@fft)
    % efficient computing for fft transform
    A = fft(A,[],3);
    if strcmp(opt,'skinny') == 1 || strcmp(opt,'econ') == 1
        min12 = min(n1,n2);
        U = zeros(n1,min12,n3);
        S = zeros(min12,min12,n3);
        V = zeros(n2,min12,n3);
        
        halfn3 = ceil((n3+1)/2);
        for i = 1 : halfn3
            [U(:,:,i),S(:,:,i),V(:,:,i)] = svd(A(:,:,i),'econ');
        end
        for i = halfn3+1 : n3
            U(:,:,i) = conj(U(:,:,n3+2-i));
            V(:,:,i) = conj(V(:,:,n3+2-i));
            S(:,:,i) = S(:,:,n3+2-i);
        end        
        if strcmp(opt,'skinny') == 1
            s1 = diag(sum(S,3))/n3^2;
            tol = max(n1,n2)*eps(max(s1));
            trank = sum(s1 > tol); % tensor tubal rank
            U = U(:,1:trank,:);
            V = V(:,1:trank,:);
            S = S(1:trank,1:trank,:);
        end
        
    elseif strcmp(opt,'full') == 1
        U = zeros(n1,n1,n3);
        S = zeros(n1,n2,n3);
        V = zeros(n2,n2,n3);
        
        halfn3 = ceil((n3+1)/2);
        for i = 1 : halfn3
            [U(:,:,i),S(:,:,i),V(:,:,i)] = svd(A(:,:,i));
        end
        for i = halfn3+1 : n3
            U(:,:,i) = conj(U(:,:,n3+2-i));
            V(:,:,i) = conj(V(:,:,n3+2-i));
            S(:,:,i) = S(:,:,n3+2-i);
        end
    end
    U = ifft(U,[],3);
    S = ifft(S,[],3);
    V = ifft(V,[],3);
else
    % other transform
    A = lineartransform(A,transform);
    if strcmp(opt,'skinny') == 1 || strcmp(opt,'econ') == 1
        min12 = min(n1,n2);
        U = zeros(n1,min12,n3);
        S = zeros(min12,min12,n3);
        V = zeros(n2,min12,n3);
        for i = 1 : n3
            [U(:,:,i),S(:,:,i),V(:,:,i)] = svd(A(:,:,i),'econ');
        end        
        if strcmp(opt,'skinny') == 1
            if transform.l > 0
                % property L'*L=L*L'=l*I, for some l>0, holds.
                s1 = diag(sum(S,3))/n3/transform.l;
            else
                % the above property does not hold for transform L
                s1 = diag(sum(S,3))/n3/norm(transform.L)^2;
            end
            tol = 1e-10;
            trank = sum(s1 > tol); % tensor tubal rank
            U = U(:,1:trank,:);
            V = V(:,1:trank,:);
            S = S(1:trank,1:trank,:);
        end
    elseif strcmp(opt,'full') == 1
        U = zeros(n1,n1,n3);
        S = zeros(n1,n2,n3);
        V = zeros(n2,n2,n3);
        for i = 1 : n3
            [U(:,:,i),S(:,:,i),V(:,:,i)] = svd(A(:,:,i));
        end
    end
    U = inverselineartransform(U,transform);
    S = inverselineartransform(S,transform);
    V = inverselineartransform(V,transform);
end

