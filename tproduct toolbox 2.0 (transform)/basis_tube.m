function ek = basis_tube(k,n3,transform)

% The tube basis under real linear transform, a tensor of size 1*1*n3
%   For fft, the tube basis ek is a tensor with the (1,1,k)-th entry
%       equaling 1 and the rest equaling 0, 
%       i.e., [ek]_(1,1,k)=1 and 0 otherwise.
%
%   For general real linear transform, the tube basis ek is a tensor
%       satisfying that L(ek) is with the (1,1,k)-th entry
%       equaling 1 and the rest equaling 0,
%       i.e., [L(ek)]_(1,1,k) = 1 and 0 otherwise.
%
% The definition of tube basis for fft can be found at our works [1][2].
% The definition of tube basis for general real linear transforms can be found at our work [3].
%
% [1] Canyi Lu, et al., Exact Low Tubal Rank Tensor Recovery from Gaussian Measurements. IJCAI. 2018
% [2] Canyi Lu, et al., Tensor Robust Principal Component Analysis with A New Tensor Nuclear Norm, TPAMI, 2019
% [3] Canyi Lu, et al., Exact Recovery of Tensor Robust Principal Component Analysis under Linear Transforms. arXiv preprint arXiv:1907.08288. 2019
%
%
% Remark: 
%       The definitions of the tube basis for fft and general real linear 
%       transforms are different. If fft is used as the linear transform L,
%       the resulted tube basis under L is different from the one defined 
%       for fft. The tube basis is not defined for other complex transforms.  
%       See more discussions in our work [3].
%
%
%
% Input:
%       k       -   tube index of nonzero entry
%       n3      -   tube basis is a tensor of size 1*1*n3
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
% Output: ek
%
%
%
% See also basis_column, unit_eijk, lineartransform, inverselineartransform
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

if k > n3
    error('%d out of bound %d.', k, n3);
end

if nargin < 3
    % fft is the default transform
    transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
end

ek = zeros(1,1,n3);
if isequal(transform.L,@fft)
    ek(1,1,k) = 1; % fft transform
elseif isa(transform.L,'function_handle')
    ek(1,1,k) = 1;
    ek = inverselineartransform(ek,transform);
elseif ismatrix(transform.L)
    if isreal(transform.L)
        % L is a real matrix
        ek(1,1,k) = 1;
        ek = inverselineartransform(ek,transform);
    else
        % L is a complex matrix
        error('the basis for complex matrix transform is not defined.');
    end
end

