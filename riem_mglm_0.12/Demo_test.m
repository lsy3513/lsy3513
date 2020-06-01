% DEMO on SPD manifolds
clear;
addpath(genpath(pwd))
rng(222);
disp('Start.')
for i = 1:10
    seed = randi(1000);
    rng(seed);
    
    Yp = randspd(3,2,3);
    V = randn(size(Yp,1));
    V = (V+V')/2*0.1;
    %Y = expmap_spd(Yp,V)
    %Y2 = addnoise_spd(Yp,0.1,10)
    error = norm_TpM_spd(Yp,V);   
    fprintf('normal %f \n', ...
        error)

    V = trnd(5,size(Yp,1),size(Yp,1));
    V = (V+V')/2;
    error = norm_TpM_spd(Yp,V);  
    fprintf('t %f \n', ...
        error)
end