function eijk = unit_eijk(i,j,k,n1,n2,n3,transform)

% Return the n1*n2*n3 sized tensor eijk which is related to the unit tensor under linear transform
%   For fft, eijk is the unit tensor with its (i,j,k)-th entry
%       equaling 1 and the rest equaling 0, 
%       i.e., [eijk]_(i,j,k)=1 and 0 otherwise.
%
%
%   For general real linear transform, eijk satisfies that L(eijk) is the 
%       unit tensor with its (i,j,k)-th entry
%       equaling 1 and the rest equaling 0, 
%       i.e., [L(eijk)]_(i,j,k)=1 and 0 otherwise.
%
% The definition and property of the unit tensor for fft can be found at our works [1][2].
% The definition and property of the unit tensor for general real linear transforms can be found at our work [3].
%
% [1] Canyi Lu, et al., Exact Low Tubal Rank Tensor Recovery from Gaussian Measurements. IJCAI. 2018
% [2] Canyi Lu, et al., Tensor Robust Principal Component Analysis with A New Tensor Nuclear Norm, TPAMI, 2019
% [3] Canyi Lu, et al., Exact Recovery of Tensor Robust Principal Component Analysis under Linear Transforms. arXiv preprint arXiv:1907.08288. 2019
%
%
% Remark:
%       Tensor eijk has a property that eijk=ei*ek*ej^T, where ei and ej 
%       are column basis and ek is the tube basis. Note that the 
%       definitions of tube basis for fft and general real linear transform
%       are different. So the definitions of eijk for fft and general real 
%       linear transforms are different. If fft is used as the linear 
%       transform L, the resulted eijk under L is different from the one 
%       defined for fft. Tensor eijk is not defined for other complex 
%       transforms. See more discussions in our work [3].
%
%
%
% Input:
%     i,j,k     -   index of nonzero entry of eijk or L(eijk)
%   n1,n2,n3    -   eijk size n1*n2*n3
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
% Output: eijk
%
%
%
% See also basis_column, basis_tube, lineartransform, inverselineartransform
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

if i > n1
    error('%d out of bound %d.', i, n1);
elseif j > n2
    error('%d out of bound %d.', j, n2);
elseif k > n3
    error('%d out of bound %d.', k, n3);
end

if nargin < 7
    % fft is the default transform
    transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
end
eijk = zeros(n1,n2,n3);
if isequal(transform.L,@fft)
    eijk(i,j,k) = 1; % fft transform
elseif isa(transform.L,'function_handle')
    eijk(i,j,k) = 1;
    eijk = inverselineartransform(eijk,transform);
elseif ismatrix(transform.L)
    if isreal(transform.L)
        % L is a real matrix
        eijk(i,j,k) = 1;
        eijk = inverselineartransform(eijk,transform);
    else
        % L is a complex matrix
        error('the basis for complex matrix transform is not defined.');
    end
end

