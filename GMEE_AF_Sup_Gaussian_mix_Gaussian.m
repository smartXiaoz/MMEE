clear all; close all;
%% MEE-AF-V10
p = 10;
q = 10;            % 滑动窗口法中的窗口长度（窗口长有利于误差曲线的平缓）
L = 4000+q;
iter = 50;         % 算法迭代次数
lamda2 = 0.996;     % RMEE
lamda3 = 0.992;     % RLS
% lamda2 = 1;       % RMEE for speech signal
lamda1 = 0.942;     % RMC

tic
for mm = 1:iter
%     vv = randn(1,L) * 0.1;
%     vv = rand(1,L) -0.5;
    v1=randn(1,L)*0.1; v2=randn(1,L)*10;
    rp=rand(1,L);
    %    vv = v1 + (rp>0.95).*v2;
    %% 设置噪声
     vv = (rp<=0.95).*v1 + (rp>0.95).*v2;   % 脉冲噪声(超高斯)
%     vv = randn(1,L) * 1;                    % 高斯噪声
%     s_sub=fun_generate_subsig5(L);
%     vv = s_sub(4,:)*5;                       % 脉冲噪声(次高斯)       
    
    
    
    %     vG = exp(-vv.^2/2/sigma1^2).*vv;
    %      vv = exp(-vv.^2/2/sigma1^2).*vv;
    
   wo1 = randn(p,1);wo2 = randn(p,1);wo3 = randn(p,1);
%       wo1 = [0,0,0.9,0,0,0,0.2,0,0,0]';
%     wo = [ kron(wo1, ones(1,L/3)) kron(wo2, ones(1,L/3)) kron(wo3, ones(1,L/3))];
     wo = [ kron(wo1, ones(1,L)) ];
    % wo = [ kron(randn(p,1), ones(1,L/2)),  kron(randn(p,1), ones(1,L/2)) ];
    uu = randn(p,L);
%     path(path,'e:\work\speech enhancement\dbs');
%     [s1,FS,NBITS]=wavread('sp02.wav');
%     s1 = s1*1;
%     
%     for ii = 1 : p
%         uu(ii, :) = s1(3000+ii : 3000+ii+L-1)*10;
%     end
    for ii = 1 : L
        dd(ii) = wo(:,ii)' * uu(:,ii) + vv(ii);
    end
    
    % dd = wo' * uu + vv;
    
    w_LMS = randn(p,1);   
    w_RLS = w_LMS;
    w_LMF = w_LMS;
    w_C_RLS =w_LMS;  %%MCC
    w_GC_LMS =w_LMS;
    w_GE_LMS =w_LMS;
%     w_M_RLS =w_LMS;%%MEE
    %% LMS
    mu1 = 0.015;
    for ii = 1 : L
        Err_LMS(mm,ii) = (wo(:,ii) - w_LMS)' *  (wo(:,ii) - w_LMS);
        dn = dd(ii);
        un = uu(:,ii);
        en = dn - w_LMS' * un;
        w_LMS = w_LMS + mu1 * (en) * un;
        %  w_LMS = w_LMS + mu1 * tanh(en) * un;
    end
        %% LMF
    mu_LMF = 0.00006;
    for ii = 1 : L
        Err_LMF(mm,ii) = (wo(:,ii) - w_LMF)' *  (wo(:,ii) - w_LMF);
        dn = dd(ii);
        un = uu(:,ii);
        en = dn - w_LMF' * un;
        w_LMF = w_LMF + mu_LMF * (en^3) * un;
        %  w_LMS = w_LMS + mu1 * tanh(en) * un;
    end
    %% RLS-MCC
    sigma1= 0.39;        % RMC
    Pn = eye(p)*1.6;
    for ii = 1 : L
        Err_MCC_RLS(mm,ii) = (wo(:,ii)  - w_C_RLS)' * (wo(:,ii)  - w_C_RLS);
        dn = dd(ii);
        un = uu(:,ii);
        en = dn - w_C_RLS' * un;
        
        kn = Pn * un / ( exp(en^2/2/sigma1^2)*lamda1 + un' * Pn * un );
        %  kn = Pn * un / ( lamda+ un' * Pn * un );
        Pn = 1/lamda1 * ( Pn - kn * un' * Pn);
        w_C_RLS = w_C_RLS +kn * en;
    end  
    %% GMCC_LMS
    alpha_c = 3;
%     lambda = 1/(beta^alpha);
    lambda_c = 0.015;
    mu_gc = 0.019;     % 步长越小，稳态误差越小，收敛越慢
    for ii = 1 : L
        Err_GMCC_LMS(mm,ii) = (wo(:,ii) - w_GC_LMS)' *  (wo(:,ii) - w_GC_LMS);
        dn = dd(ii);
        un = uu(:,ii);
        en = dn - w_GC_LMS' * un;
        exp_mid = exp(-1*lambda_c*(abs(en))^alpha_c);
        w_GC_LMS = w_GC_LMS + mu_gc*exp_mid*(abs(en))^(alpha_c-1)*sign(en)*un;
%         w_GMCC_LMS = w_GMCC_LMS + mu1 * (en) * un;
        %    w_LMS = w_LMS + mu1 * tanh(en) * un;
    end
    
     
%% GMEE_LMS
    mu_ge = 0.075;
    alpha_e = 1;                    % 越小精度越高，收敛越慢
    beta_e = 100;
    lambda_e =beta_e^alpha_e; 
%     lambda_e =8;                   % 越大收敛越快，貌似不影响精度           
    
    for ii = q : L
        Err_GMEE_LMS(mm,ii) = (wo(:,ii) - w_GE_LMS)' *  (wo(:,ii) - w_GE_LMS);
        % 计算该窗口内的所有样本的误差
        for jj = 1 : q
            error_q(jj) = dd(ii - q + jj) - w_GE_LMS' * uu(: , ii - q + jj);
        end
        u_q = uu(:,ii - q + 1:ii);     % 一个窗口内的输入
        sum_sum = 0;
        for i = 1:q
            for j=1:q
                G = (alpha_e/2*beta_e*(gamma(1/alpha_e)))*exp(-1*(abs(error_q(i)-error_q(j)))^alpha_e/lambda_e);
                ei_ej = (abs(error_q(i)-error_q(j)))^(alpha_e-1);
                sum_sum = sum_sum + G*ei_ej*(u_q(:,i)-u_q(:,j))*sign(error_q(i)-error_q(j));
            end
        end
        w_GE_LMS = w_GE_LMS + mu_ge*(alpha_e/(q^2*beta_e^alpha_e))*sum_sum;
    end
end
toc
figure(1),hold on;
plot(10* log10(mean(Err_LMS(:,1:L-q))),'g','LineWidth',1.5),
plot(10* log10(mean(Err_LMF(:,1:L-q))),'LineWidth',1.5),
plot(10* log10(mean(Err_GMCC_LMS(:,1:L-q))),'r','LineWidth',1.5),
plot(10* log10(mean(Err_GMEE_LMS(:,q:L))),'c','LineWidth',1.5),
% plot(10* log10(mean(Err_MCC_RLS(:,1:L-q))),'LineWidth',1.5),
legend('LMS','LMF','GMCC','GMEE');
% legend('LMS','LMF','GMCC','RMC');
axis([0,L-q,-35,20]);
xlabel('Iterations');ylabel('MSD');

E_LMS = 10* log10(mean(Err_LMS));
E_LMF = 10* log10(mean(Err_LMF));
E_GMCC_LMS = 10* log10(mean(Err_GMCC_LMS));
E_GMEE_LMS = 10* log10(mean(Err_GMEE_LMS));
E_MCC_RLS = 10* log10(mean(Err_MCC_RLS));

E_LMS = mean(E_LMS(end/2:end));
E_LMF = mean(E_LMF(end/2:end));
E_GMCC_LMS = mean(E_GMCC_LMS(end/2:end));
E_GMEE_LMS = mean(E_GMEE_LMS(end/2:end));
E_MCC_RLS = mean(E_MCC_RLS(end/2:end));


disp(['E_LMS:',num2str(E_LMS)]);
disp(['E_LMF:',num2str(E_LMF)]);
disp(['E_GMCC_LMS:',num2str(E_GMCC_LMS)]);
disp(['E_GMEE_LMS:',num2str(E_GMEE_LMS)]);
disp(['E_MCC_RLS:',num2str(E_MCC_RLS)]);


