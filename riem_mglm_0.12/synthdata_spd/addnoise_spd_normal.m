% 점 p 에 대해 err가 isotropic하게 새로운 점 q를 반환

function q = addnoise_spd_normal(p, err_scale, maxerr)
    n = size(p,1);
    err_size = randn_f(n,err_scale,maxerr);
    x = randpoint(n);
    expD = diag(exp(err_size * x));
    U = randortho(n);
    expS = U * expD * U';
    rtp = sqrtm(p);
    q = rtp * expS * rtp;
end

% S = randS(size(p,1));
% rtp = sqrtm(p);
% q = rtp * expm(err_scale * S) * rtp;

% n 차원 unit sphere 에서 uniform random하게 한 점 뽑기
function p = randpoint(n)
    x = randn(n,1);
    p = x/sqrt(sum(x .* x));
end

% uniform random 하게 n * n orthonormal matrix 뽑기
function ortho = randortho(n)
    X = zeros(n);
    for i = 1 : n
        X(:,i) = randpoint(n);
    end
    ortho = gram(X);
end

% uniform random 하게 tr(S^2) = 1 인 S 뽑기
function S = randS(n)
    lambda = randpoint(n);
    D = diag(lambda);
    U = randortho(n);
    S = U * D * U';
end

% 그람-슈미트 방법. 인수는 n * n 행렬, 반환인자도 n * n 행렬이고, 각 열 벡터에 대해 그람슈미트 방법을 적용.
function Y = gram(X)
    Y = zeros(size(X));
    for k = 1 : size(X,1)
        proj_temp = zeros(size(X,1),1);
        for j = 1:(k-1)
            proj_temp = proj_temp + proj_uv(Y(:,j),X(:,k));
        end
        Y(:,k) = X(:,k) - proj_temp;
    end
    for k = 1 : size(X,1)
        Y(:,k) = Y(:,k)/norm(Y(:,k));
    end
end

function proj = proj_uv(u,v)
    proj = sum(u .* v)/sum(u .* u) * u;
end

% sample from e(-x^2/2c^2) * x^(n-1) (0 < x < maxerr)
%(n=dimY=3)
function x = randn_f(n,c,maxerr) 
    while true
        while true
            x = abs(randn)*c;
            if x < maxerr
                break
            end
        end
        if rand < (x/maxerr)^(n-1)
            break
        end
    end
end