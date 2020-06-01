% addpath(genpath(pwd))
% 
for i = 1:size(Y_temp,4)
    Y = Y_temp(:,:,:,i);
    Ybar = karcher_mean_spd(Y,[],500);
    [pls, Vls, Els, Yhatls, gnormls] = mglm_spd(X_temp,Y);
    r2_iterative  = r2stat_spd(Ybar, Y, Yhatls);
    fprintf("%f\n",r2_iterative)
end

% mse_ls = mse_for_RealData_spd_ls(X_temp,Y_temp(:,:,:,1));
% mse_huber = mse_for_RealData_spd_huber(X_temp,Y_temp(:,:,:,1));
% mse_l1 = mse_for_RealData_spd_l1(X_temp,Y_temp(:,:,:,1));