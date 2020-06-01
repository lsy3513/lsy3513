% DEMO on the unit sphere.
clear;
addpath(genpath(pwd));
rng(22)
disp('Start.');

for iter = 1:10

    mse_p_ls = 0;
    mse_p_hb = 0;
    mse_p_l1 = 0;
    mse_p_tk = 0;

    mse_v_ls = 0;
    mse_v_hb = 0;
    mse_v_l1 = 0;
    mse_v_tk = 0;

    seed = randi(1000);
    %rng(seed);
    %fprintf('seed: %f\n', seed)

    synth_sphere_data_mse

    for i = 1:mse_iter
        Y = Y_raw(:,:,i);
        Ybar = karcher_mean_sphere(Y,[],500);   

        [pls, Vls, Els, Yhatls, gnormls,gradp_list] = mglm_sphere(X,Y);
        mse_p_ls = mse_p_ls + norm(logmap_sphere(Yp(:,1),pls))^2;
        Vls_par = paralleltranslateAtoB_sphere(pls, Yp(:,1),Vls);
        for j = 1:npivots
            mse_v_ls = mse_v_ls + norm(V(:,j)-Vls_par(:,j))^2;
        end

        [phb, Vhb, Ehb, Yhathb, gnormhb] = mglm_sphere_huber(X,Y);
        mse_p_hb = mse_p_hb + norm(logmap_sphere(Yp(:,1),phb))^2;
        Vhb_par = paralleltranslateAtoB_sphere(phb, Yp(:,1),Vhb);
        for j = 1:npivots
            mse_v_hb = mse_v_hb + norm(V(:,j)-Vhb_par(:,j))^2;
        end
        
        [pl1, Vl1, El1, Yhatl1, gnorml1] = mglm_sphere_l1(X,Y);
        mse_p_l1 = mse_p_l1 + norm(logmap_sphere(Yp(:,1),pl1))^2;
        Vl1_par = paralleltranslateAtoB_sphere(pl1, Yp(:,1),Vl1);
        for j = 1:npivots
            mse_v_l1 = mse_v_l1 + norm(V(:,j)-Vl1_par(:,j))^2;
        end
% 
%         [ptk, Vtk, Etk, Yhattk, gnormtk] = mglm_sphere_tukey(X,Y);
%         mse_p_tk = mse_p_tk + norm(logmap_sphere(Yp(:,1),ptk))^2;
%         Vtk_par = paralleltranslateAtoB_sphere(ptk, Yp(:,1),Vtk);
%         for j = 1:npivots
%             mse_v_tk = mse_v_tk + norm(V(:,j)-Vtk_par(:,j))^2;
%         end

        %x = [Y(1,:),Yhatls(1,:),Yhathb(1,:),Yhattk(1,:)];
        %y = [Y(2,:),Yhatls(2,:),Yhathb(2,:),Yhattk(2,:)];
        %labels = cellstr( num2str([1:400]') );  %' # labels correspond to their order

        %plot(x, y, 'rx')
        %text(x, y, labels, 'VerticalAlignment','bottom','HorizontalAlignment','right')
    end
    %fprintf('mse_p_ls %f mse_p_hb %f mse_p_l1 %f mse_p_tk %f\n', mse_p_ls/mse_iter, mse_p_hb/mse_iter, mse_p_l1/mse_iter, mse_p_tk/mse_iter)
    %fprintf('mse_v_ls %f mse_v_hb %f mse_V_l1 %f mse_v_tk %f\n', mse_v_ls/mse_iter, mse_v_hb/mse_iter, mse_v_l1/mse_iter, mse_v_tk/mse_iter)
    fprintf('%f,%f,%f,%f,', mse_p_ls/mse_iter, mse_p_hb/mse_iter, mse_p_l1/mse_iter, mse_p_tk/mse_iter)
    fprintf('%f,%f,%f,%f\n', mse_v_ls/mse_iter, mse_v_hb/mse_iter, mse_v_l1/mse_iter, mse_v_tk/mse_iter)
end