function ei = basis_column(i,n,n3,transform)

% The column basis under real linear transform, a tensor of size n*1*n3
%   For fft, the column basis ei is a tensor with the (i,1,1)-th entry
%       equaling 1 and the rest equaling 0, 
%       i.e., [ei]_(i,1,1)=1 and 0 otherwise.
%
%   For general real linear transform, the column basis ei is a tensor
%       satisfying that L(ei) is with the (i,1,:) tube being an all one
%       tube and the rest equaling 0, 
%       i.e., [L(ei)]_(i,1,:) = 1 and 0 otherwise.
%
% The definition of column basis for fft can be found at our works [1][2].
% The definition of column basis for general real linear transforms can be found at our work [3].
%
% [1] Canyi Lu, et al., Exact Low Tubal Rank Tensor Recovery from Gaussian Measurements. IJCAI. 2018
% [2] Canyi Lu, et al., Tensor Robust Principal Component Analysis with A New Tensor Nuclear Norm, TPAMI, 2019
% [3] Canyi Lu, et al., Exact Recovery of Tensor Robust Principal Component Analysis under Linear Transforms. arXiv preprint arXiv:1907.08288. 2019
%
%
% Remark: 
%       The definition of the column basis for general real linear 
%       transforms is a generalization of the column basis for fft. If fft is 
%       used as the linear transform L, the resulted column basis under L 
%       reduces to the one defined for fft. The column basis is not defined 
%       for other complex transforms. See more discussions in our work [3].
%
%
%
% Input:
%       i       -   row index of nonzero entries
%       n, n3   -   column basis is a tensor of size n*1*n3
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
% Output: ei
%
%
%
% See also basis_tube, unit_eijk, lineartransform, inverselineartransform
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

if i > n
    error('%d out of bound %d.', i, n);
end

if nargin < 4
    % fft is the default transform
    transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
end

ei = zeros(n,1,n3);
if isequal(transform.L,@fft)
    ei(i,1,1) = 1; % fft transform
elseif isa(transform.L,'function_handle')
    ei(i,1,:) = ones(1,n3);
    ei = inverselineartransform(ei,transform);
elseif ismatrix(transform.L)
    if isreal(transform.L)
        % L is a real matrix
        ei(i,1,:) = ones(1,n3);
        ei = inverselineartransform(ei,transform);
    else
        % L is a complex matrix
        error('the basis for complex matrix transform is not defined.');
    end
end
