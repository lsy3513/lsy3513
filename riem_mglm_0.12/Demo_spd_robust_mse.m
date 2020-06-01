% DEMO on SPD manifolds
clear;
addpath(genpath(pwd))
%rng(2);
disp('Start.')

for iter = 1:5

    mse_p_ls = 0;
    mse_p_hb = 0;
    mse_p_l1 = 0;
    mse_p_tk = 0;

    mse_v_ls = 0;
    mse_v_hb = 0;
    mse_v_l1 = 0;
    mse_v_tk = 0;

    seed = randi(1000);

    synth_dti_data_mse

    for i = 1:mse_iter
        Y = Y_raw(:,:,:,i);
        Ybar = karcher_mean_spd(Y,[],500);   

        [pls, Vls, Els, Yhatls, gnormls] = mglm_spd(X,Y);
        mse_p_ls = mse_p_ls + norm_TpM_spd(Yp(:,:,1),logmap_spd(Yp(:,:,1),pls))^2;
        Vls_par = paralleltranslateAtoB_spd(pls, Yp(:,:,1),Vls);
        for j = 1:npivots
            mse_v_ls = mse_v_ls + norm_TpM_spd(Yp(:,:,1),V(:,:,j)-Vls_par(:,:,j))^2;
        end

        [phb, Vhb, Ehb, Yhathb, gnormhb] = mglm_spd_huber(X,Y);
        mse_p_hb = mse_p_hb + norm_TpM_spd(Yp(:,:,1),logmap_spd(Yp(:,:,1),phb))^2;
        Vhb_par = paralleltranslateAtoB_spd(phb, Yp(:,:,1),Vhb);
        for j = 1:npivots
            mse_v_hb = mse_v_hb + norm_TpM_spd(Yp(:,:,1),V(:,:,j)-Vhb_par(:,:,j))^2;
        end
         
        [pl1, Vl1, El1, Yhatl1, gnorml1] = mglm_spd_l1(X,Y);
        mse_p_l1 = mse_p_l1 + norm_TpM_spd(Yp(:,:,1),logmap_spd(Yp(:,:,1),pl1))^2;
        Vl1_par = paralleltranslateAtoB_spd(pl1, Yp(:,:,1),Vl1);
        for j = 1:npivots
            mse_v_l1 = mse_v_l1 + norm_TpM_spd(Yp(:,:,1),V(:,:,j)-Vl1_par(:,:,j))^2;
        end

        % [ptk, Vtk, Etk, Yhattk, gnormtk] = mglm_spd_tukey(X,Y);
        % mse_p_tk = mse_p_tk + norm_TpM_spd(Yp(:,:,1),logmap_spd(Yp(:,:,1),ptk))^2;
        % Vtk_par = paralleltranslateAtoB_spd(ptk, Yp(:,:,1),Vtk);
        % for j = 1:npivots
        %     mse_v_tk = mse_v_tk + norm_TpM_spd(Yp(:,:,1),V(:,:,j)-Vtk_par(:,:,j))^2;
        % end

    end
    fprintf('%f,%f,%f,%f,', mse_p_ls/mse_iter, mse_p_hb/mse_iter, mse_p_l1/mse_iter, mse_p_tk/mse_iter)
    fprintf('%f,%f,%f,%f\n', mse_v_ls/mse_iter, mse_v_hb/mse_iter, mse_v_l1/mse_iter, mse_v_tk/mse_iter)
end
