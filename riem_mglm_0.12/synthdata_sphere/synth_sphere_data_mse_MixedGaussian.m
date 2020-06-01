%% Parameters
dimY = 15;
noise = 0.8;
npairs = 10;
udist = 3; % udist : V의 크기를 결정하는 parameter
noutliers = 1;
ntest = 10;
npivots = 5;

%% For figure
%X = [0:0.25:1 0:0.25:1; 0 0 0 0 0 1 1 1 1 1  ];
X = rand(npivots,npairs);
X = center(X);

%% Synthesized data
npivots = size(X,1); % Number of points except the base point
Yp = zeros(dimY, npivots + 1);
while true
    Yp(:,1) = unitvec(randn(dimY,1));
    Yp(:,2) = addnoise_sphere(Yp(:,1),1.5,0.5);
    Yp(:,3) = addnoise_sphere(Yp(:,1),1.5,0.5);
    V = zeros(dimY,npivots);
    for j = 1:npivots
        V(:,j) = logmap_sphere(Yp(:,1),Yp(:,j+1));
    end
    
    %% Generate Ground Truth
    Y_0 = zeros(dimY, size(X,2));
    issafe = true;
    for i = 1:length(X)
        Vtmp = V*X(:,i);
        if norm(Vtmp) > udist
            issafe = false;
            break
        end
        Y_0(:,i) = expmap_sphere(Yp(:,1),Vtmp);
    end

    if issafe
        break
    end
end

%% Add noise and make some samples.
Xsample = [];
Ysample = zeros(dimY,size(Y_0,2)*npairs,ntest);

for k = 1:ntest
    isample = 1;
    for i = 1:(npairs - noutliers)
        for j = 1:size(Y_0,2)
            Ysample(:,isample,k) = addnoise_sphere(Y_0(:,j),noise,0.017);
            isample = isample + 1;
        end
    end
    for i = (npairs - noutliers+1):npairs
        for j = 1:size(Y_0,2)
            Ysample(:,isample,k) = addnoise_sphere(Y_0(:,j),noise,0.274);
            isample = isample + 1;
        end
    end
end

for i = 1:npairs
    Xsample =[Xsample X];
end

X = Xsample;
Y_raw = unitvec(Ysample);
