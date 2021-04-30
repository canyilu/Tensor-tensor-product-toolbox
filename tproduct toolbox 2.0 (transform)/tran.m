function At = tran(A,transform)

% The conjugate transpose of a 3 way tensor under linear transform
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
%
% Output: At    -   n2*n1*n3 tensor
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
At = zeros(n2,n1,n3);
if isequal(transform.L,@fft)
    % fft transform
    At(:,:,1) = A(:,:,1)';
    for i = 2 : n3
        At(:,:,i) = A(:,:,n3-i+2)';
    end       
elseif isa(transform.L,'function_handle')
    A = lineartransform(A,transform);
    for i = 1 : n3
        At(:,:,i) = A(:,:,i)';
    end
    At = inverselineartransform(At,transform);
elseif ismatrix(transform.L)
    if isreal(transform.L)
        % L is a real matrix
        for i = 1 : n3
            At(:,:,i) = A(:,:,i)';
        end
    else
        % L is a complex matrix
        A = lineartransform(A,transform);
        for i = 1 : n3
            At(:,:,i) = A(:,:,i)';
        end
        At = inverselineartransform(At,transform);
    end
end

