A = [1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
n = size(A,1);
err_scale=0.3;
maxerr = 1.0;
for iter = 1:2000
    r2 = rand_norm_old(A,err_scale,maxerr);
    r1 = sphere_randn(n,err_scale,maxerr);
    fprintf("%f,%f\n",r1,r2);
end
% f(k) ~ exp(-k^2/2c^2) * sin(k) (0 < k < maxerr < pi) and 0 (o.w) 를 만족하게 x 뽑음
function x = sphere_randn(n,c,maxerr)
    while true
        while true
            k = abs(randn * c);
            if k < maxerr
                break
            end
        end
        if rand < sin(k)^(n-2)
            x = k;
            break
        end
    end
end


function r = rand_norm_old(A,c,maxerr)
    V = randn(size(A))*c;
    V = V-V'*A*A;
    r = norm(V);
    if norm(V) > maxerr
        r = maxerr;
    end
end