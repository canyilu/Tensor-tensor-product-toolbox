function [Q,R] = tqr(A,transform,opt)

% Tensor orthogonal-triangular decomposition under linear transform
%
%   [Q,R] = tqr(A,transform), where A is n1*n2*n3, produces an n1*n2*n3 upper triangular
%   tensor R and an n1*n1*n3 orthogonal tensor Q so that A = Q*R.
%
%   [Q,R] = tqr(A,transform,'econ') produces the "economy size" decomposition.
%   If n1>n2, only the first n2 lateral slices of Q and the first n2 horizontal
%   slices of R are computed. If n1<=n2, this is the same as [Q,R] = tqr(A).
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
%      opt      -   options for different outputs of Q and R:
%                  - 'econ': produces the "economy size" decomposition. 
%                  - If not specified (default), produces full decomposition, i.e., A = Q*R, where
%                     A - n1*n2*n3
%                     Q - n1*n1*n3
%                     R - n1*n2*n3
%
%
% Output: Q, R
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

[n1,n2,n3] = size(A);
if nargin < 2
    % fft is the default transform
    transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
end

if isequal(transform.L,@fft)
    % efficient computing for fft transform
    A = fft(A,[],3);    
    if n1>n2 && exist('opt', 'var') && strcmp(opt,'econ') == 1
        Q = zeros(n1,n2,n3);
        R = zeros(n2,n2,n3);
        halfn3 = ceil((n3+1)/2);
        for i = 1 : halfn3
            [Q(:,:,i),R(:,:,i)] = qr(A(:,:,i),0);
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
            [Q(:,:,i),R(:,:,i)] = qr(A(:,:,i));
        end
        for i = halfn3+1 : n3
            Q(:,:,i) = conj(Q(:,:,n3+2-i));
            R(:,:,i) = conj(R(:,:,n3+2-i));
        end
    end
    Q = ifft(Q,[],3);
    R = ifft(R,[],3);
else
    % other transform
    A = lineartransform(A,transform);
    if n1>n2 && exist('opt', 'var') && strcmp(opt,'econ') == 1
        Q = zeros(n1,n2,n3);
        R = zeros(n2,n2,n3);
        for i = 1 : n3
            [Q(:,:,i),R(:,:,i)] = qr(A(:,:,i),0);
        end
    else
        Q = zeros(n1,n1,n3);
        R = zeros(n1,n2,n3);
        for i = 1 : n3
            [Q(:,:,i),R(:,:,i)] = qr(A(:,:,i));
        end
    end
    Q = inverselineartransform(Q,transform);
    R = inverselineartransform(R,transform);
end
