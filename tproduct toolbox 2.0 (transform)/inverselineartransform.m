function A = inverselineartransform(A,transform)

% perform inverse transform of the given transform on tensor A along 3rd dim
%
%   A           - n1*n2*n3 tensor
%   transform   - a structure which defines the linear transform
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
%
%
%	Examples: 
%            n = 10; A = rand(n,n,n); A = inverselineartransform(A,transform);
%            Here transform is defined below:
%              
%           - linear transform as the funtion handle:
%             (1) fft: transform.L = @fft, transform.l = n3, transform.inverseL = @ifft. 
%                      A = inverselineartransform(A,transform) is equivalent to A = ifft(A,[],3);
%             (2) dct: transform.L = @dct, transform.l = 1, transform.inverseL = @idct. 
%                      A = inverselineartransform(A,transform) is equivalent to A = idct(A,[],3);
%           - linear transform as the invertible matrix: transform.L = L, transform.l = l, transform.inverseL = inverseL, where
%             (1) matrix corresponding to discrete Fourier transform (equivalent to fft): L = dftmtx(n3); l = n3; inverseL = L'/l;
%             (2) matrix corresponding to discrete cosine transform (equivalent to dct):  L = dct(eye(n3)); l = 1; inverseL = L';
%             (3) random orthogonal real matrix: L = RandOrthMat(n3); l = 1; inverseL = L';
%             (4) any invertible matrix, i.e., L = randn(n3); l = -1; inverseL = inv(L);
%
%
%
% See also lineartransform
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

if nargin < 2
    A = ifft(A,[],3); % fft is the default transform
elseif isa(transform.L,'function_handle')
    A = transform.inverseL(A,[],3);
elseif ismatrix(transform.L)
    n3 = size(A,3);
    [l1,l2] = size(transform.L);
    if l1 ~= l2 || l1 ~= n3
        error('Inner tensor dimensions must agree.');
    end
    if transform.l > 0
        % property L'*L=L*L'=l*I, for some l>0, holds. Then, inv(L)=L'/l
        A = tmprod(A,transform.L'/transform.l,3);
    else
        % arbitrary invertible matrix without the above property
        if isfield(transform,'inverseL')
            A = tmprod(A,transform.inverseL,3);
        else
            A = tmprod(A,inv(transform.L),3);
        end
    end
end

