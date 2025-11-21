clc;
clear; close all;
addpath('Test_data/','Functions/');
addpath(genpath('Tucker/'));
rand('seed',60);
randn('seed',60);

% Image selection
fig={'baboon'};
% Sampling ratio set
SR_set=[0.1 0.2 0.3];

for i=1:3
    fig1=fig{1};
    imgfile1 = strcat('Test_data/',fig1,'.bmp');
    X = double(imread(imgfile1));
    X1 = X;
    dim = size(X1);
    % Generate known data
    SR = SR_set(i);
    P = round(SR*prod(dim));
    known = randsample(prod(dim),P);
    [known,~] = sort(known);
    % Observation
    data = X1(known);
    
    %Tucker
    %TMac-rank_inc strategy
    disp('begining TMac');
    opts.alpha_adj = 0;
    opts.rank_adj = 1*ones(1,3);
    opts.rank_min = [1 1 1];
    opts.rank_max = [100,100,3];
    EstCoreNway = [5 5 1];
    Ntucker = 3;
    coNway = zeros(1,Ntucker);
    Nwaytucker = [256,256,3];
    for n = 1:Ntucker
        coNway(n) = prod(Nwaytucker)/Nwaytucker(n);
    end
    for m = 1:3
        X0{m} = randn(Nwaytucker(m),EstCoreNway(m));
        Y0{m} = randn(EstCoreNway(m),coNway(m));
    end
    opts.X0 = X0; opts.Y0 = Y0;
    [X_inc,Y_inc,Out_inc] = TMac(data,known,Nwaytucker,EstCoreNway,opts);       
    Mrec = zeros(Nwaytucker);
    for m = 1:Ntucker
        Mrec = Mrec+Out_inc.alpha(m)*Fold(X_inc{m}*Y_inc{m},Nwaytucker,m);
    end
    SSIM_TMac = ssim(Mrec,X1);
    PSNR_TMac=psnr_stdc(X1,Mrec);
    disp(['Image_',fig1,'_SR_',num2str(SR), '_SSIM_',num2str(SSIM_TMac),'_PSNR_', num2str(PSNR_TMac)])

    % TMac-ARE
    disp('begining TMac-ARE');
    [TMacARE_Result] = Tucker_ARE(X1,X1,known);
    SSIM_TMacARE = ssim(TMacARE_Result.Xc,X1);
    PSNR_TMacARE = psnr_stdc(X1,TMacARE_Result.Xc);
    disp(['Image_',fig1,'_SR_',num2str(SR), '_SSIM_',num2str(SSIM_TMacARE),'_PSNR_', num2str(PSNR_TMacARE)])
end

