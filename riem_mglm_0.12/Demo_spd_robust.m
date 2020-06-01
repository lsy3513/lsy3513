% DEMO on SPD manifolds
clear;
addpath(genpath(pwd))
rng(223);
disp('Start.')

mse_p_ls = 0;
mse_p_hb = 0;
mse_p_tk = 0;

mse_v_ls = 0;
mse_v_hb = 0;
mse_v_tk = 0;

for i = 1:3
    seed = randi(1000);
    rng(seed);
    
    synth_dti_data3

    Ybar = karcher_mean_spd(Y,[],500);

    [pls, Vls, Els, Yhatls, gnormls] = mglm_spd(X,Y);
    r2_iterative_ls  = r2stat_spd(Ybar, Y, Yhatls);
    mse_p_ls = mse_p_ls + norm_TpM_spd(Yp(:,:,1),logmap_spd(Yp(:,:,1),pls))^2;
    for j = 1:size(V,3)
        mse_v_ls = mse_v_ls + norm_TpM_spd(Yp(:,:,1),V(:,:,j)-Vls(:,:,j))^2;
    end

    [phb, Vhb, Ehb, Yhathb, gnormhb] = mglm_spd_huber(X,Y);
    r2_iterative_hb  = r2stat_spd(Ybar, Y, Yhathb);
    mse_p_hb = mse_p_hb + norm_TpM_spd(Yp(:,:,1),logmap_spd(Yp(:,:,1),phb))^2;
    for j = 1:size(V,3)
        mse_v_hb = mse_v_hb + norm_TpM_spd(Yp(:,:,1),V(:,:,j)-Vhb(:,:,j))^2;
    end

    [ptk, Vtk, Etk, Yhattk, gnormtk] = mglm_spd_tukey(X,Y);
    r2_iterative_tk  = r2stat_spd(Ybar, Y, Yhattk);
    mse_p_tk = mse_p_tk + norm_TpM_spd(Yp(:,:,1),logmap_spd(Yp(:,:,1),ptk))^2;
    for j = 1:size(V,3)
        mse_v_tk = mse_v_tk + norm_TpM_spd(Yp(:,:,1),V(:,:,j)-Vtk(:,:,j))^2;
    end

    fprintf('seed %-4d: r2_iterative_ls %f r2_iterative_hb %f r2_iterative_tk %f\n', seed, r2_iterative_ls, r2_iterative_hb, r2_iterative_tk)
end
fprintf('mse_p %f mse_p_hb %f mse_p_tk %f\n', mse_p_ls, mse_p_hb, mse_p_tk)
fprintf('mse_v %f mse_v_hb %f mse_v_tk %f\n', mse_v_ls, mse_v_hb, mse_v_tk)
    %[ple, Vle, Ele, Yhatle, logY, Yv ] = mglm_logeuc_spd(X,Y);

    %r2_logeuc  = r2stat_spd(Ybar, Y, Yhatle);
    %fprintf('seed %-4d: r2_iterative %f r2_logeuc %f\n', seed, r2_iterative, r2_logeuc)