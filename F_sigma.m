function [steay_mee,steay_mee_hampel,steay_mee_andrew,steay_mee_huber] = F_sigma(sigma)
Epoch = 12;
LEN = 10000;
L = 15;
m = 5;
mu_mee = 0.001;
mu_andrews = 0.001;
mu_huber = 0.001;
mu_hampel = 0.001;
err_mee = zeros(Epoch, LEN-L+1);
err_hampel = zeros(Epoch, LEN-L+1);
err_andrew = zeros(Epoch, LEN-L+1);
err_huber = zeros(Epoch, LEN-L+1);

parfor mm = 1 : Epoch
    e_greedy = rand(1, LEN);%均匀分布
    wo = randn(m, 1);%初始w0为30*1
    UU = randn(m, LEN);%uu为输入，理论上为30*1，但取3000组数据；
    DD = zeros(LEN,1);
    %% case1
    vv=randn(1,LEN)*0.1;
    v1 = randn(1,LEN) * 0.1;
    v2 = randn(1,LEN) * 10;
    VV = (e_greedy <= 0.95).*v1 + (e_greedy >0.95).*v2; %non-gaussian-noise      
    %% case5 mutimodel
    v1 = randn(1,LEN) * 0.3 + 0.2;
    v2 = randn(1,LEN) * 20;
    v3 = randn(1,LEN) * 0.3 -0.2;
    VV = (e_greedy <= 0.4).*v1 + ((e_greedy >0.4) & (e_greedy <0.6)).*v2 + (e_greedy >= 0.6).*v3; %non-gaussian-noise
    for ii = 1 : LEN
        DD(ii) = wo' * UU(:,ii) + VV(ii);%产生1*3000的行向量，也即是输出
%         DD(ii) = wo' * UU(:,ii) +alpha_noise(ii);%产生1*3000的行向量，也即是输出
    end
    w_lms = randn(m, 1);    
    err_mee(mm,:) = F_MEE(LEN,mu_mee,wo,w_lms,DD,UU,sigma,L);
    err_hampel(mm,:) = F_MEE_hampel(LEN,mu_hampel,wo,w_lms,DD,UU,sigma,L);
    err_andrew(mm,:) = F_MEE_andrew(LEN,mu_andrews,wo,w_lms,DD,UU,sigma,L);
    err_huber(mm,:) = F_MEE_huber(LEN,mu_huber,wo,w_lms,DD,UU,sigma,L);
end

steay_mee = mean(mean(err_mee(:,end-3000:end)));
steay_mee_hampel = mean(mean(err_hampel(:,end-3000:end)));
steay_mee_andrew = mean(mean(err_andrew(:,end-3000:end)));
steay_mee_huber = mean(mean(err_huber(:,end-3000:end)));