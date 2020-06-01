function Anew = addnoise_sphere_t(A, nu, err_scale, maxerr)
    % 길이 dim(Y)(=15) 짜리 unit vector를 받아
        n = size(A,1);
        V = randn(size(A)); %randn(n,m) 은 표준정규분포 난수 n*m 행렬, size(A)는 (15,1) 반환, 즉 여기서 V 는 N(0,1)로 이루어진 15 * 1 벡터.
        V = V-V'*A*A; % T_A(M)으로 V를 projection, 즉 V에서 A와 수직인 성분 분해한 거임.
        V = sphere_randt(n,nu, err_scale,maxerr) * V / norm(V);
        Anew = expmap_sphere(A, V);
    end

    % f(k) ~ trnd(k;nu,c) * sin(k)^(n-2) (0 < k < maxerr < pi) and 0 (o.w) 를 만족하게 x 뽑음
    function x = sphere_randt(n,nu,c,maxerr)
        while true
            while true
                x = abs(trnd(nu))*c;
                if x < maxerr
                    break
                end
            end
            if maxerr >= pi/2
                M = 1;
            else
                M = sin(maxerr)^(n-2);
            end
            if rand < sin(x)^(n-2)/M
                break
            end
        end
    end
    

    % 밑은 과거의 유물
    % f(k) ~ (1+k^2/(c^2*nu))^(-(nu+1)/2) * sin(k)^(n-2) (0 < k < maxerr < pi) and 0 (o.w) 를 만족하게 x 뽑음
    % 이는 f(point)~ (1+k^2/nu)^(-(nu+1)/2) 가 되기 위함임.
    % k = sqrt((n-2)*nu/(nu-n+3)) 에서 최대. nu > n-3이 아니면, maxerr가 최대
    % 주의 n>=3 가정임.
    % function x = sphere_randt(n,nu,c,maxerr)
    %     if nu > n-3
    %         k_max = sqrt((n-2)*nu/(nu-n+3));
    %         M = (1+k_max^2/nu)^(-(nu+1)/2) * k_max^(n-2);
    %     end
    %     if nu <= n-3
    %         if maxerr >= 1
    %             M = (1+1/nu)^(-(nu+1)/2);
    %         end
    %         if maxerr < 1
    %             k_max = maxerr;
    %             M = (1+k_max^2/nu)^(-(nu+1)/2) * k_max^(n-2);
    %         end
    %     end
    %     while true
    %         k = rand * maxerr;
    %         if rand < (1+k^2/nu)^(-(nu+1)/2) * sin(k)^(n-2) / M
    %             x = k;
    %             break
    %         end
    %     end
    % end
    
    % 밑은 구 버전
    % function x = sphere_randn(n,c,maxerr)
    %     while true
    %         while true
    %             k = abs(randn * c);
    %             if k < maxerr
    %                 break
    %             end
    %         end
    %         if rand < sin(k)^(n-2)
    %             x = k;
    %             break
    %         end
    %     end
    % end
    
        
    