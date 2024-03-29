%mee -lms
clear all;
close all;
LEN = 2000;
m = 5;
L = 15;

mu_lms = 0.02;
mu_lmm = 0.02;
mu_gmcc = 0.019;

mu_mee = 0.0001;
mu_mee0 = 0.0022;
mu_andrews = 0.0002;
mu_huber = 0.0002;
mu_hampel = 0.0003;
% mu_andrews = 0.001;
% mu_huber = 0.001;
% mu_hampel = 0.001;
%case3
% mu_lms = 0.02;
% mu_lmm = 0.02;
% mu_gmcc = 0.019;
% mu_mee = 0.4;
% mu_andrews = 0.3;
% mu_huber = 0.3;
% mu_hampel = 0.3;

sigma0 = 0.8; 
sigma1 = 4; 
sigma2 = 5; 
% sigma2 = 0.8; 
tic
Epoch = 20;
err_lms0 = zeros(Epoch, LEN);
err_lms11 = zeros(Epoch, LEN);
Err_GMCC_LMS = zeros(Epoch, LEN);
err_lms = zeros(Epoch, LEN-L+1);
err_lms1 = zeros(Epoch, LEN-L+1);
err_lms6 = zeros(Epoch, LEN-L+1);
err_lms5 = zeros(Epoch, LEN-L+1);
parfor mm = 1 :Epoch
% err_lms0 = zeros(Epoch, LEN);
% err_lms11 = zeros(Epoch, LEN);
% Err_GMCC_LMS = zeros(Epoch, LEN);
% err_lms = zeros(Epoch, LEN);
% err_lms1 = zeros(Epoch, LEN);
% err_lms6 = zeros(Epoch, LEN);
% err_lms5 = zeros(Epoch, LEN);
    e_greedy = rand(1, LEN);
    %vv = randn(1,LEN) * 0.1;
    wo = randn(m, 1);
    UU = randn(m, LEN);
    DD = zeros(LEN,1);
    %% case1
    vv=randn(1,LEN)*0.1;
    v1 = randn(1,LEN) * 0.1;
    v2 = randn(1,LEN) * 10;
    VV = (e_greedy <= 0.95).*v1 + (e_greedy >0.95).*v2; %non-gaussian-noise
    
    %% case2

    alpha_noise = alpha_stable_noise(1.2,0.1,0,0,LEN);%alpha=1.2 gamma=0.1
    varNoise = 1;
    VV = sqrt(varNoise) * alpha_noise;

    %% case3 chi-r
    VV = chi2rnd(1,[1, LEN]) * 0.1;


    %% case4 t_dis
    VV = trnd(1.5,[1, LEN]) * 0.1;
   %% case5 mutimodel
%     v1 = randn(1,LEN) * 0.3 + 0.2;
%     v2 = randn(1,LEN) * 20;
%     v3 = randn(1,LEN) * 0.3 -0.2;
%     VV = (e_greedy <= 0.4).*v1 + ((e_greedy >0.4) & (e_greedy <0.6)).*v2 + (e_greedy >= 0.6).*v3; %non-gaussian-noise
    % case6 mutimodel
%     v1 = randn(1,LEN) * 0.3 - 0.2;
%     v2 = randn(1,LEN) * 10;
%     v3 = trnd(1.5,[1, LEN]) * 0.1 + 0.2;
%     VV = (e_greedy <= 0.8).*v1 + ((e_greedy >0.8) & (e_greedy <0.9)).*v2 + (e_greedy >= 0.9).*v3; %non-gaussian-noise
%     
    for ii = 1 : LEN
        DD(ii) = wo' * UU(:,ii) + VV(ii);
    end

    w_lms = randn(m, 1);
    
    en = zeros(L,1);%L=10��
    
     
    Err_GMCC_LMS(mm,:) = F_GCCLMS(LEN,mu_gmcc,wo,w_lms,DD,UU);
    err_lms0(mm,:) = F_LMS(LEN,mu_lms,wo,w_lms,DD,UU);
    err_lms11(mm,:) = F_LMM(LEN,mu_lmm,wo,w_lms,DD,UU);  
    err_lms(mm,:) = F_MEE(LEN,mu_mee,wo,w_lms,DD,UU,sigma1,L);
    err_mee0(mm,:) = F_MEE(LEN,mu_mee0,wo,w_lms,DD,UU,sigma0,L);
    err_lms5(mm,:) = F_MEE_hampel(LEN,mu_hampel,wo,w_lms,DD,UU,sigma2,L);
    err_lms1(mm,:) = F_MEE_andrew(LEN,mu_andrews,wo,w_lms,DD,UU,sigma2,L);
    err_lms6(mm,:) = F_MEE_huber(LEN,mu_huber,wo,w_lms,DD,UU,sigma2,L);
 
   
    
    
    disp(mm)
end
toc


figure(1),hold on;
box on;
wid = 2;
MarkerSize = 8;
Indices = 100;


color_lms0 = [0.85 0.33 0.10]; 
color_lms11 = [0.93 0.69 0.13]; 
color_gmcc = [0.49 0.18 0.56]; 
color_lms = [0.6350 0.0780 0.1840]; 

plot(10*log10(mean(err_lms0)), '-', 'color', color_lms0, 'LineWidth', wid, 'Marker', 'o', 'MarkerSize', MarkerSize, 'MarkerFaceColor', color_lms0, 'MarkerIndices', 1:Indices:LEN)
plot(10*log10(mean(err_lms11)), '-', 'color', color_lms11, 'LineWidth', wid, 'Marker', 's', 'MarkerSize', MarkerSize, 'MarkerFaceColor', color_lms11, 'MarkerIndices', 1:Indices:LEN)
plot(10*log10(mean(Err_GMCC_LMS)), '-', 'color', color_gmcc, 'LineWidth', wid, 'Marker', '^', 'MarkerSize', MarkerSize, 'MarkerFaceColor', color_gmcc, 'MarkerIndices', 1:Indices:LEN)
plot(10*log10(mean(err_lms)), '-c', 'LineWidth', wid, 'Marker', 'd', 'MarkerSize',MarkerSize , 'MarkerFaceColor', 'c', 'MarkerIndices', 1:Indices:LEN)
plot(10*log10(mean(err_mee0)), '-m', 'LineWidth', wid, 'Marker', 'd', 'MarkerSize',MarkerSize , 'MarkerFaceColor', 'm', 'MarkerIndices', 1:Indices:LEN)

plot(10*log10(mean(err_lms1)), '-r', 'LineWidth', wid, 'Marker', 'p', 'MarkerSize',MarkerSize , 'MarkerFaceColor', 'r', 'MarkerIndices', 1:Indices:LEN)
plot(10*log10(mean(err_lms6)), '-g', 'LineWidth', wid, 'Marker', 'h', 'MarkerSize', MarkerSize, 'MarkerFaceColor', 'g', 'MarkerIndices', 1:Indices:LEN)
plot(10*log10(mean(err_lms5)), '-b', 'LineWidth', wid, 'Marker', '*', 'MarkerSize', MarkerSize, 'MarkerFaceColor', 'b', 'MarkerIndices', 1:Indices:LEN)


h=legend(['LMS ({\eta\rm}=' num2str(mu_lms) ') '],['LMM (\eta=' num2str(mu_lmm) ', \Delta_1=3, \Delta_2=4, \Delta_3=5)']...
,['GMCC ({\eta=' num2str(mu_gmcc) '}, \alpha=3, \lambda=0.015)'],['MEE0 ({\eta=' num2str(mu_mee) '}, L=15, \sigma=' num2str(sigma1) ')'],['MEE1 ({\eta=' num2str(mu_mee0) '}, L=15, \sigma=' num2str(sigma0) ')']...
,['MMEE-Andrews ({\eta=' num2str(mu_andrews) '}, L=15, \sigma=' num2str(sigma2) ', \Delta_1=1)'],['MMEE-Huber ({\eta=' num2str(mu_huber) '}, L=15, \sigma=' num2str(sigma2) ', \Delta_1=0.75)']...
,['MMEE-Hample ({\eta=' num2str(mu_hampel) '}, L=15, \sigma=' num2str(sigma2) ', \Delta_1=0.5, \Delta_2=2, \Delta_3=4)']);
set(h,'FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',24)
xlabel('Iteration');
ylabel('MSD');


