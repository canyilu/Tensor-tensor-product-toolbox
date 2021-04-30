function tnn = tnn(A,transform)

% Tensor nuclear norm of a 3 way tensor under linear transform
%
% An error message is printed if the property L'*L=L*L'=l*I, for some l>0,
% of the linear transform does not hold.
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
% Output: 
%        tnn    -   tensor nuclear norm
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

n3 = size(A,3);
if nargin == 2
    if transform.l < 0
        error("property L'*L=L*L'=l*I does not holds for some l>0.");
    end
else    
    % fft is the default transform
    transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
end

tnn = 0;
if isequal(transform.L,@fft)
    % efficient computing for fft transform
    A = fft(A,[],3);
    % i=1
    s = svd(A(:,:,1),'econ');
    tnn = tnn+sum(s);    
    % i=2,...,halfn3
    halfn3 = round(n3/2);
    for i = 2 : halfn3
        s = svd(A(:,:,i),'econ');
        tnn = tnn+sum(s)*2;
    end    
    % if n3 is even
    if mod(n3,2) == 0
        i = halfn3+1;
        s = svd(A(:,:,i),'econ');
        tnn = tnn+sum(s);
    end
    tnn = tnn/n3;
else
    A = lineartransform(A,transform);
    for i = 1 : size(A,3)
        s = svd(A(:,:,i),'econ');
        tnn = tnn+sum(s);
    end
    tnn = tnn/transform.l;
end
