function Anew = addnoise_sphere_old(A,c,maxerr)
    V = randn(size(A))*c;
    V = V-V'*A*A;
    if norm(V) > maxerr
        V = V/norm(V)*maxerr;
    end
    Anew = expmap_sphere(A, V);
    