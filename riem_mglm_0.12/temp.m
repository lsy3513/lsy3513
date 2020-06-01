clc
clear
fi0_plt = 0;
fi1_plt = 360;
R_plt = 1;
R0_plt = 0;
R1_plt = 1;
M_plt = 30;
dfi_plt = (fi1_plt-fi0_plt)/M_plt;
dR_plt = (R1_plt-R0_plt)/M_plt;
fi_plt = [fi0_plt:dfi_plt:fi1_plt];
aa_plt=pi/180;
theta0_plt=0;
theta1_plt=360;
dtheta_plt=(theta1_plt-theta0_plt)/M_plt;
theta_plt=[theta0_plt:dtheta_plt:theta1_plt];
t_plt=0;
tt_plt=-1;
for j=1:M_plt/2+1
    tt_plt=tt_plt+1;
    t_plt=tt_plt;
    for i=1:M_plt+1
        t_plt=t_plt+1;
        b_plt(t_plt)=t_plt;
        x_plt(i)=R_plt*sin(aa_plt*fi_plt(i))*cos(aa_plt*theta_plt(j));
        y_plt(i)=R_plt*sin(aa_plt*fi_plt(i))*sin(aa_plt*theta_plt(j));
        z_plt(i)=R_plt*cos(aa_plt*fi_plt(i));
        %pause
        %v(1:6,b(t))=[ x(i) y(i) z(i) fi(i) theta(j) b(t)];
    end
plot3(x_plt,y_plt,z_plt)
axis equal
grid on
hold on
end

t_plt=0;
tt_plt=-1;
for j=1:M_plt/2+1
    tt_plt=tt_plt+1;
    t_plt=tt_plt;
    for i=1:M_plt+1
        t_plt=t_plt+1;
        bb(t_plt)=t_plt;
        xx_plt(i)=R_plt*sin(aa_plt*fi_plt(j))*cos(aa_plt*theta_plt(i));
        yy_plt(i)=R_plt*sin(aa_plt*fi_plt(j))*sin(aa_plt*theta_plt(i));
        zz_plt(i)=R_plt*cos(aa_plt*fi_plt(j));
        %pause
        %ang=45;
        %RR = rotx(ang);
        %vv = [xx(i);yy(i);zz(i)];
        %pause
        %yyy = RR*vv;
        %vvv(1:6,bb(t))=[ yyy(1) yyy(2) yyy(3) fi(j) theta(i) bb(t)];
        plot3(xx_plt,yy_plt,zz_plt)
        %plot3(vvv(1,bb(t)),vvv(2,bb(t)),vvv(3,bb(t)),'r*')
        axis equal
        grid on
        hold on
    end
end

addpath(genpath(pwd));
rng(222)
disp('Start.');

for iter = 1:3

    mse_p_ls = 0;
    mse_p_hb = 0;
    mse_p_tk = 0;

    mse_v_ls = 0;
    mse_v_hb = 0;
    mse_v_tk = 0;

    seed = randi(1000);
    %rng(seed);
    fprintf('seed: %f\n', seed)

    synth_sphere_data_mse
    
    plot3(Yp(1,1),Yp(2,1),Yp(3,1),'r*')

    for i = 1:ntest
        Y = Y_raw(:,:,i);
        Ybar = karcher_mean_sphere(Y,[],500);   

        [pls, Vls, Els, Yhatls, gnormls] = mglm_sphere(X,Y);
        mse_p_ls = mse_p_ls + norm(logmap_sphere(Yp(:,1),pls))^2;
        for j = 1:npivots
            mse_v_ls = mse_v_ls + norm(V(:,j)-Vls(:,j))^2;
        end

        [phb, Vhb, Ehb, Yhathb, gnormhb] = mglm_sphere_huber(X,Y);
        mse_p_hb = mse_p_hb + norm(logmap_sphere(Yp(:,1),phb))^2;
        for j = 1:npivots
            mse_v_hb = mse_v_hb + norm(V(:,j)-Vhb(:,j))^2;
        end

        [ptk, Vtk, Etk, Yhattk, gnormtk] = mglm_sphere_tukey(X,Y);
        mse_p_tk = mse_p_tk + norm(logmap_sphere(Yp(:,1),ptk))^2;
        for j = 1:npivots
            mse_v_tk = mse_v_tk + norm(V(:,j)-Vtk(:,j))^2;
        end
        
        plot333(pls,phb,ptk)

        %x = [Y(1,:),Yhatls(1,:),Yhathb(1,:),Yhattk(1,:)];
        %y = [Y(2,:),Yhatls(2,:),Yhathb(2,:),Yhattk(2,:)];
        %labels = cellstr( num2str([1:400]') );  %' # labels correspond to their order

        %plot(x, y, 'rx')
        %text(x, y, labels, 'VerticalAlignment','bottom','HorizontalAlignment','right')
    end
    fprintf('mse_p_ls %f mse_p_hb %f mse_p_tk %f\n', mse_p_ls/ntest, mse_p_hb/ntest, mse_p_tk/ntest)
    fprintf('mse_v_ls %f mse_v_hb %f mse_v_tk %f\n', mse_v_ls/ntest, mse_v_hb/ntest, mse_v_tk/ntest)
end

function plot333(point1,point2,point3)
    plot3(point1(1,1),point1(2,1),point1(3,1),'b*')
    plot3(point2(1,1),point2(2,1),point2(3,1),'g*')
    plot3(point3(1,1),point3(2,1),point3(3,1),'m*')
end