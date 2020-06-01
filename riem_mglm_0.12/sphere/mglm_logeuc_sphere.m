function [p, V, E, Y_hat, U] = mglm_logeuc_sphere(X,Y,varargin)
%MGLM_SPHERE performs MGLM on the unit sphere by by Log Euclidean
%framework.
%
%   [p, V, E, Y_hat, U] = MGLM_LOGEUC_SPHERE(X, Y)
%
%   [p, V, E, Y_hat, gnorm] = MGLM_SPHERE(X, Y, niter)
%   has optional parameter niter for Karcher mean calculation.  
%
%   The result is in p, V, E, Y_hat.
%
%   X is dimX x N column vectors
%   Y is dimY x N column vectors (points on the unit sphere in R^dimY).
%
%   p is a base point.
%   V is a set of tangent vectors (dimY x dimX).
%   E is the sum of squared geodesic error.
%   Y_hat is the prediction.
%   U is the orthogonal basis of TpM.
%
%   See also MGLM_SPHERE, EXPMAP_SPHERE, FEVAL_SPHERE, PREDICTION_SPHERE,
%   LOGMAP_VECS_SPHERE, PARALLELTRANSLATEATOB_SPHERE

%   Hyunwoo J. Kim
%   $Revision: 0.1 $  $Date: 2014/06/23 17:50:37 $

if nargin >=3
    niter = varargin{1};
else
    niter = 500;
end
% Linear transform

[ ndim ndata] = size(X);
Xc = X - repmat(mean(X,2),1,ndata); % mean(X,2) 는 각 행벡터의 평균값을 반환한 열벡터, Xc는 중앙화 시킨 값
p = karcher_mean_sphere(Y, ones(ndata,1)/ndata, niter);
logY = logmap_vecs_sphere(p, Y);

% Get orthogonal bases of TpM
U = null(ones(size(p,1),1)*p'); % null(A) 는 A의 null space의 orthonormal basis 를 열벡터로 갖는 행렬


Yu = U'*logY; %logY is represented by U
% Yu = L*X
L = Yu/Xc; % LXc = Yu, 이거 항상 있음?
V = U*L;
logY_hat = V*Xc;
Y_hat = expmap_vecs_sphere(p,logY_hat);
E = gsqerr_sphere(Y, Y_hat); % Geodesic squared error