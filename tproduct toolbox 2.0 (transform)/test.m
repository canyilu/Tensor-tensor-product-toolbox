%
% References:
%
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
% Written by Canyi Lu (canyilu@gmail.com)
%

clear
%% tprod
n1 = 2;
n2 = 3;
n3 = 5;
m2 = 3;
A = rand(n1,n2,n3);
B = rand(n2,m2,n3);

% transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
% transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
% transform.L = dftmtx(n3); transform.l = n3;
% transform.L = dct(eye(n3)); transform.l = 1;
% transform.L = RandOrthMat(n3); transform.l = 1;
transform.L = randn(n3); transform.l = -1; transform.inverseL = inv(transform.L); % arbitrary invertible matrix
transform.L = randn(n3)+i*(rand(n3)-1); transform.l = -1; transform.inverseL = inv(transform.L); % arbitrary invertible matrix

C = tprod(A,B,transform);

%% tran
n1 = 2;
n2 = 2;
n3 = 2;
A = rand(n1,n2,n3);

% transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
% transform.L = dftmtx(n3); transform.l = n3;
% transform.L = dct(eye(n3)); transform.l = 1;
% transform.L = RandOrthMat(n3); transform.l = 1;
% transform.L = randn(n3); transform.l = -1; transform.inverseL = inv(transform.L); % arbitrary invertible matrix
transform.L = randn(n3)+i*rand(n3); transform.l = -1; transform.inverseL = inv(transform.L); % arbitrary invertible matrix

At = tran(A,transform);

%% teye
n1 = 3;
n2 = 3;
n3 = 4;
A = rand(n1,n2,n3);

% transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
% transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
% transform.L = dftmtx(n3); transform.l = n3;
% transform.L = dct(eye(n3)); transform.l = 1;
% transform.L = RandOrthMat(n3); transform.l = 1;
% transform.L = randn(n3); transform.l = -1; transform.inverseL = inv(transform.L); % arbitrary invertible matrix
transform.L = randn(n3)+i*rand(n3); transform.l = -1; transform.inverseL = inv(transform.L); % arbitrary invertible matrix

I = teye(n1,n3,transform);
AI1 = tprod(A,I,transform);
AI2 = tprod(I,A,transform);
dif_teye1 = norm(A(:)-AI1(:))
dif_teye2 = norm(A(:)-AI2(:))

%% tubalrank
n1 = 15;
n2 = 15;
n3 = 200;
r = 10;
A = rand(n1,r,n3);
B = rand(r,n2,n3);

% transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
% transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
% transform.L = dftmtx(n3); transform.l = n3;
% transform.L = dct(eye(n3)); transform.l = 1;
% transform.L = RandOrthMat(n3); transform.l = 1;
transform.L = randn(n3); transform.l = -1; transform.inverseL = inv(transform.L); % arbitrary invertible matrix
transform.L = randn(n3)+i*rand(n3); transform.l = -1; transform.inverseL = inv(transform.L); % arbitrary invertible matrix


C = tprod(A,B,transform);
trank = tubalrank(C,transform);
dif_tubalrank = trank - r

%% tsn
n1 = 30;
n2 = 40;
n3 = 20;
A = rand(n1,n2,n3);

% transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
% transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
% transform.L = dftmtx(n3); transform.l = n3;
% transform.L = dct(eye(n3)); transform.l = 1;
% transform.L = RandOrthMat(n3); transform.l = 1;
transform.L = randn(n3); transform.l = -1; transform.inverseL = inv(transform.L); % arbitrary invertible matrix
transform.L = randn(n3)+i*rand(n3); transform.l = -1; transform.inverseL = inv(transform.L); % arbitrary invertible matrix

[U,S,V] = tsvd(A,transform,'skinny');

tsn1 = tsn(A,transform);
Sbar = lineartransform(S,transform);
tsn2 = max(Sbar(:))
dif_tsn = tsn1 - tsn2

%% tnn
n1 = 30;
n2 = 40;
n3 = 20;
A = rand(n1,n2,n3);

% transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
% transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
% transform.L = dftmtx(n3); transform.l = n3;
transform.L = dct(eye(n3)); transform.l = 1;

[U,S,V] = tsvd(A,transform,'skinny');

tnn1 = tnn(A,transform);
tnn2 = S.*teye(n1,n3,transform); tnn2 = sum(tnn2(:));
dif_tnn = tnn1-tnn2

%% prox_tnn
n1 = 5;
n2 = 6;
n3 = 4;
A = rand(n1,n2,n3);

% transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
% transform.L = dftmtx(n3); transform.l = n3;
% transform.L = dct(eye(n3)); transform.l = 1;
% transform.L = RandOrthMat(n3); transform.l = 1;

rho = 0.5;
prox = prox_tnn(A,rho,transform);

%% tinv
n1 = 9
n2 = 9
n3 = 10
A = rand(n1,n2,n3);

% transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
% transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
% transform.L = dftmtx(n3); transform.l = n3;
% transform.L = dct(eye(n3)); transform.l = 1;
% transform.L = RandOrthMat(n3); transform.l = 1;
transform.L = randn(n3); transform.l = -1; transform.inverseL = inv(transform.L); % arbitrary invertible matrix
transform.L = randn(n3)+i*rand(n3); transform.l = -1; transform.inverseL = inv(transform.L); % arbitrary invertible matrix

Ainv = tinv(A,transform);
Ieye = teye(n1,n3,transform);
AAinv1 = tprod(A,Ainv,transform);
AAinv2 = tprod(Ainv,A,transform);

res_inv1 = norm(Ieye(:)-AAinv1(:))
res_inv2 = norm(Ieye(:)-AAinv2(:))

%% tsvd
n1 = 20;
n2 = 27;
n3 = 40;
r = 10;

% transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
% transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
% transform.L = dftmtx(n3); transform.l = n3;
% transform.L = dct(eye(n3)); transform.l = 1;
% transform.L = RandOrthMat(n3); transform.l = 1;
transform.L = randn(n3); transform.l = -1; transform.inverseL = inv(transform.L); % arbitrary invertible matrix
transform.L = randn(n3)+i*rand(n3); transform.l = -1; transform.inverseL = inv(transform.L); % arbitrary invertible matrix

A = rand(n1,r,n3);
B = rand(r,n2,n3);
C = tprod(A,B,transform);

[U,S,V] = tsvd(C,transform,'skinny');
USV1 = tprod(U,S,transform);
USV1 = tprod(USV1,tran(V,transform),transform);
dif_tsvd1 = norm(C(:)-USV1(:))
dif_trank = r - size(U,2)

[U,S,V] = tsvd(C,transform,'econ');
USV2 = tprod(U,S,transform);
USV2 = tprod(USV2,tran(V,transform),transform);
dif_tsvd2 = norm(C(:)-USV2(:))

[U,S,V] = tsvd(C,transform,'full');
USV3 = tprod(U,S,transform);
USV3 = tprod(USV3,tran(V,transform),transform);
dif_tsvd3 = norm(C(:)-USV3(:))

%% qr

n1 = 20;
n2 = 15;
n3 = 12;
A = rand(n1,n2,n3);

% transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
% transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
% transform.L = dftmtx(n3); transform.l = n3;
% transform.L = dct(eye(n3)); transform.l = 1;
% transform.L = RandOrthMat(n3); transform.l = 1;
% transform.L = randn(n3); transform.l = -1; transform.inverseL = inv(transform.L); % arbitrary invertible matrix
transform.L = randn(n3)+i*rand(n3); transform.l = -1; transform.inverseL = inv(transform.L); % arbitrary invertible matrix


[Q,R] = tqr(A,transform,'econ');
QR = tprod(Q,R,transform);
dif_qr1 = norm(A(:)-QR(:))

[Q,R] = tqr(A,transform);
QR = tprod(Q,R,transform);
dif_qr2 = norm(A(:)-QR(:))

%% basis_tube, basis_column, unit_eijk

n1 = 5;
n2 = 4;
n3 = 4;

i = 2;
j = 2;
k  = 2;

% transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
% transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
% transform.L = dftmtx(n3); transform.l = n3;
% transform.L = dct(eye(n3)); transform.l = 1;
% transform.L = RandOrthMat(n3); transform.l = 1;
transform.L = randn(n3); transform.l = -1; transform.inverseL = inv(transform.L); % arbitrary invertible matrix


ek = basis_tube(k,n3,transform);
Lek = lineartransform(ek,transform);

ei = basis_column(i,n1,n3,transform);
Lei = lineartransform(ei,transform);
ej = basis_column(j,n2,n3,transform);

eijk = tprod(ei,ek,transform);
eijk = tprod(eijk,tran(ej,transform),transform);
Leijk = lineartransform(eijk,transform);

eijk2 = unit_eijk(i,j,k,n1,n2,n3,transform);
Leijk2 = lineartransform(eijk2,transform);



