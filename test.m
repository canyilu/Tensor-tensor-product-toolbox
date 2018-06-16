
% References: 
% Canyi Lu, Tensor-Tensor Product Toolbox. Carnegie Mellon University. 
% June, 2018. https://github.com/canyilu/tproduct.
%
% Canyi Lu, Jiashi Feng, Yudong Chen, Wei Liu, Zhouchen Lin and Shuicheng
% Yan, Tensor Robust Principal Component Analysis with A New Tensor Nuclear
% Norm, arXiv preprint arXiv:1804.03728, 2018
%
% Written by Canyi Lu (canyilu@gmail.com)
%

clear

%% discrete Fouriter transformation and FFT

% FFT on vectors
n = 5;
v = rand(n,1);
vbar1 = fft(v);
F = dftmtx(n);
vbar2 = F*v;
difv = norm(vbar1-vbar2)

circv = gallery('circul',v)'; % circulant matrix
invF = F'/n; % inverse of F
difv2 = norm(F*circv*invF -diag(vbar1),'fro') % diagonal property

% FFT on tensors
n1 = 3;
n2 = 4;
n3 = 5;
A = rand(n1,n2,n3);
F = dftmtx(n3);
FI1 = kron(F,eye(n1));
FI2 = kron(inv(F),eye(n2));
diagAbar1 = FI1*bcirc(A)*FI2; % diagonal property
Abar = fft(A,[],3);
diagAbar2 = bdiag(Abar);
difA = norm(diagAbar1(:)-diagAbar2(:))


%% tprod: tensor-tensor product
n1 = 4;
n2 = 5;
n3 = 6;
r = 2;
A = rand(n1,r,n3);
B = rand(r,n2,n3);
C = tprod(A,B);

I = teye(n1,n3);
AI = tprod(I,A);
difAI = norm(A(:)-AI(:))


%% tsvd: tensor singular value decomposition
n1 = 4;
n2 = 5;
n3 = 6;
r = 2;
A = rand(n1,r,n3);
B = rand(r,n2,n3);
C = tprod(A,B);
[U1,S1,V1] = tsvd(C); % full t-SVD
C1 = tprod(tprod(U1,S1),tran(V1)); % C = U*S*V^*
difC1 = norm(C(:)-C1(:))
difU1 = tprod(tran(U1),U1) - teye(n1,n3); % U^*U = I
difU1 = norm(difU1(:))
difV1 = tprod(tran(V1),V1) - teye(n2,n3); % V^*V = I
difV1 = norm(difV1(:))

[U2,S2,V2] = tsvd(C,'econ'); % "economy size" t-SVD
C2 = tprod(tprod(U2,S2),tran(V2)); % C = U*S*V^*
difC2 = norm(C(:)-C2(:))
difU2 = tprod(tran(U2),U2) - teye(min(n1,n2),n3); % U^*U = I
difU2 = norm(difU2(:))
difV2 = tprod(tran(V2),V2) - teye(min(n1,n2),n3); % V^*V = I
difV2 = norm(difV2(:))


[U3,S3,V3] = tsvd(C,'skinny'); % skinny t-SVD
C3 = tprod(tprod(U3,S3),tran(V3)); % C = U*S*V^*
difC3 = norm(C(:)-C3(:))
r = tubalrank(C);
difU3 = tprod(tran(U3),U3) - teye(r,n3); % U^*U = I
difU3 = norm(difU3(:))
difV3 = tprod(tran(V3),V3) - teye(r,n3); % V^*V = I
difV3 = norm(difV3(:))


%% trank: tensor tubal rank
n1 = 4;
n2 = 5;
n3 = 6;
r = 2;
A = rand(n1,r,n3);
B = rand(r,n2,n3);
C = tprod(A,B);
trank = tubalrank(C);
diftrank = trank - r


%% tnn: tensor nuclear norm
n1 = 4;
n2 = 5;
n3 = 6;
A = rand(n1,n2,n3);

tnnA1 = tnn(A);

bcircA = bcirc(A);
tnnA2 = sum(svd(bcircA,'econ'))/n3;

Abar = fft(A,[],3);
bdiagAbar = bdiag(Abar);
tnnA3 = sum(svd(bdiagAbar,'econ'))/n3;

diftnn1 = tnnA2-tnnA1
diftnn2 = tnnA3-tnnA1

%% tsn: tensor spectral norm
n1 = 4;
n2 = 5;
n3 = 6;
A = rand(n1,n2,n3);
tsnA1 = tsn(A);

bcircA = bcirc(A);
tsnA2 = norm(bcircA,2);

Abar = fft(A,[],3);
bdiagAbar = bdiag(Abar);
tsnA3 = norm(bdiagAbar,2);

diftsn1 = tsnA2-tsnA1
diftsn2 = tsnA3-tsnA1


%% prox_tnn: proximal operator of tensor nuclear norm
n1 = 4;
n2 = 5;
n3 = 6;
A = rand(n1,n2,n3);

tau = 10;
X = prox_tnn(A,tau);
tubalrank(X)


%% tinv: inverse of a tensor
n = 4;
n3 = 6;
A = rand(n,n,n3);
Ainv = tinv(A);
AAinv1 = tprod(A,Ainv);
AAinv2 = tprod(Ainv,A);
I = teye(n,n3);
difinv1 = norm(AAinv1(:)-I(:))
difinv2 = norm(AAinv2(:)-I(:))


%% tqr: tensor QR factorization
n1 = 6;
n2 = 5;
n3 = 6;
A = rand(n1,n2,n3);
if n1 > n2
    [Q,R] = tqr(A,'econ');
    QtQ = tprod(tran(Q),Q);
    I = teye(n2,n3);
    diftqr1 = norm(QtQ(:)-I(:))
end

[Q,R] = tqr(A);
QtQ = tprod(tran(Q),Q);
I = teye(n1,n3);
diftqr2 = norm(QtQ(:)-I(:))

