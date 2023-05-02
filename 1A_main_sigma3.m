%mee -lms
clear all;
close all;
LEN = 1000;
m = 5;
L=15;
sigma1 = 5;
sigma2 = 10;
sigma3 = 20;
sigma4 = 100;
mu1 = 1.1;
mu2 = 3;
mu3 = 11;
mu4 = 600;
mu5 = 2;
mu6 = 9;
mu7 = 35;
mu8 = 1000;
tic
for mm = 1 : 500
    e_greedy = rand(1, LEN);%均匀分布
    %vv = randn(1,LEN) * 0.1;
    vv=randn(1,LEN)*0.1;
    v1 = randn(1,LEN) * 0.1;
    v2 = randn(1,LEN) * 10;
    VV = (e_greedy <= 0.9).*v1 + (e_greedy >0.9).*v2; %non-gaussian-noise
    
    wo = randn(m, 1);%初始w0为30*1
    UU = randn(m, LEN);%uu为输入，理论上为30*1，但取3000组数据；
    alpha_noise = alpha_stable_noise(1.1,0.1,0,0,LEN);%alpha=1.1 gamma=0.1
    varNoise = 1;
    alpha_noise = sqrt(varNoise) * alpha_noise;
    for ii = 1 : LEN
        DD(ii) = wo' * UU(:,ii) + VV(ii);%产生1*3000的行向量，也即是输出
%         DD(ii) = wo' * UU(:,ii) +alpha_noise(ii);%产生1*3000的行向量，也即是输出
    end

    w_lms = randn(m, 1);
    w_mee1 = w_lms;
    w_mee2 = w_lms;
    w_mee3 = w_lms;
    w_mee4 = w_lms;
    w_hampel1 = w_lms;
    w_hampel2 = w_lms;
    w_hampel3 = w_lms;
    w_hampel4 = w_lms;
    
    %%MEE
    
    en = zeros(L,1);%L=10；   
    for ii = L : LEN
        err_mee1(mm,ii-L + 1) = (w_mee1-wo)'*(w_mee1-wo);%mm为最外层循环，从1-200；最终就要画出此图像；
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_mee1' * UU(:,ii - L + jj);
            %en这个列向量始终存储10个误差1-10，2-11，3-12......
        end 
        for i = 1 : L
            for j = 1 : L
                uij = UU(:,ii - L + i) - UU(:,ii - L + j);
                eij = en(i) - en(j);
                tmp = 1/(L^2*sigma1^2) * exp(-(eij^2)/2/sigma1^2) * eij * uij;
                w_mee1 = w_mee1 + mu1 * tmp;%mul=0.2;sigma=2;
            end           
        end
    end
    
    %%MEE

    en = zeros(L,1);%L=10；   
    for ii = L : LEN
        err_mee2(mm,ii-L + 1) = (w_mee2-wo)'*(w_mee2-wo);%mm为最外层循环，从1-200；最终就要画出此图像；
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_mee2' * UU(:,ii - L + jj);
            %en这个列向量始终存储10个误差1-10，2-11，3-12......
        end 
        for i = 1 : L
            for j = 1 : L
                uij = UU(:,ii - L + i) - UU(:,ii - L + j);
                eij = en(i) - en(j);
                tmp = 1/(L^2*sigma2^2) * exp(-(eij^2)/2/sigma2^2) * eij * uij;
                w_mee2 = w_mee2 + mu2 * tmp;%mul=0.2;sigma=2;
            end           
        end
    end
    
    %%MEE

    en = zeros(L,1);%L=10；   
    for ii = L : LEN
        err_mee3(mm,ii-L + 1) = (w_mee3-wo)'*(w_mee3-wo);%mm为最外层循环，从1-200；最终就要画出此图像；
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_mee3' * UU(:,ii - L + jj);
            %en这个列向量始终存储10个误差1-10，2-11，3-12......
        end 
        for i = 1 : L
            for j = 1 : L
                uij = UU(:,ii - L + i) - UU(:,ii - L + j);
                eij = en(i) - en(j);
                tmp = 1/(L^2*sigma3^2) * exp(-(eij^2)/2/sigma3^2) * eij * uij;
                w_mee3 = w_mee3 + mu3 * tmp;%mul=0.2;sigma=2;
            end           
        end
    end
    
    %%MEE

    en = zeros(L,1);%L=10；   
    for ii = L : LEN
        err_mee4(mm,ii-L + 1) = (w_mee4-wo)'*(w_mee4-wo);%mm为最外层循环，从1-200；最终就要画出此图像；
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_mee4' * UU(:,ii - L + jj);
            %en这个列向量始终存储10个误差1-10，2-11，3-12......
        end 
        for i = 1 : L
            for j = 1 : L
                uij = UU(:,ii - L + i) - UU(:,ii - L + j);
                eij = en(i) - en(j);
                tmp = 1/(L^2*sigma4^2) * exp(-(eij^2)/2/sigma4^2) * eij * uij;
                w_mee4 = w_mee4 + mu4 * tmp;%mul=0.2;sigma=2;
            end           
        end
    end
    
%% robust5:Hampel 法

    en = zeros(L,1);
    for ii = L : LEN
        err_hampel1(mm,ii-L + 1) = (w_hampel1-wo)'*(w_hampel1-wo);
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_hampel1' * UU(:,ii - L + jj);
        end 
        u = 0.6745 * en / median(abs(en - median(en)));
        %median(A)是求矩阵A每列的中位数，若个数为偶数，则为最中间两个值的平均值；
        %en - median(en)：en为10*1的列向量，该向量每个元素都减去其中位数,加绝对值会使每个元素都变为非负数；
        %u的分母最终是一个数，u为en的倍数（10*1列向量）
        aa = 0.5;
        bb = 2;
        cc = 4;
        w5 = zeros(L);
        jw5 = zeros(1,L);
        for jj = 1 : L
            if abs(u(jj)) <= aa
                w5(jj) = 1;
                jw5(:,jj) = en(jj);
            elseif (aa < abs(u(jj))) && (abs(u(jj)) <=bb)
                w5(jj) = aa / abs(u(jj));
                jw5(:,jj) = en(jj) * w5(jj);
            elseif (bb < abs(u(jj))) && (abs(u(jj)) <=cc)
                w5(jj) = (aa / abs(u(jj))) * ((cc-abs(u(jj)))/(cc-bb));
                jw5(:,jj) = en(jj) * w5(jj);
            elseif abs(u(jj)) > cc
                w5(jj) = 0;
                jw5(:,jj) = 0;
            end
        end
        for i = 1 : L
            for j = 1 : L   
                %这个权值是乘在了en上，不加权重eij=en(i) - en(j);
                eij = jw5(:,i) - jw5(:,j);
                jw = UU(:,ii - L + i) - UU(:,ii - L + j);%jw就是uij；
                tmp =1/(L^2*sigma1^2) * exp(-(eij^2)/2/sigma1^2) * eij * jw;
                w_hampel1 = w_hampel1 + mu5 * tmp;
            end
        end
    end
    
    
    
    %% robust5:Hampel 法

    en = zeros(L,1);
    for ii = L : LEN
        err_hampel2(mm,ii-L + 1) = (w_hampel2-wo)'*(w_hampel2-wo);
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_hampel2' * UU(:,ii - L + jj);
        end 
        u = 0.6745 * en / median(abs(en - median(en)));
        %median(A)是求矩阵A每列的中位数，若个数为偶数，则为最中间两个值的平均值；
        %en - median(en)：en为10*1的列向量，该向量每个元素都减去其中位数,加绝对值会使每个元素都变为非负数；
        %u的分母最终是一个数，u为en的倍数（10*1列向量）
        aa = 0.5;
        bb = 2;
        cc = 4;
        w5 = zeros(L);
        jw5 = zeros(1,L);
        for jj = 1 : L
            if abs(u(jj)) <= aa
                w5(jj) = 1;
                jw5(:,jj) = en(jj);
            elseif (aa < abs(u(jj))) && (abs(u(jj)) <=bb)
                w5(jj) = aa / abs(u(jj));
                jw5(:,jj) = en(jj) * w5(jj);
            elseif (bb < abs(u(jj))) && (abs(u(jj)) <=cc)
                w5(jj) = (aa / abs(u(jj))) * ((cc-abs(u(jj)))/(cc-bb));
                jw5(:,jj) = en(jj) * w5(jj);
            elseif abs(u(jj)) > cc
                w5(jj) = 0;
                jw5(:,jj) = 0;
            end
        end
        for i = 1 : L
            for j = 1 : L   
                %这个权值是乘在了en上，不加权重eij=en(i) - en(j);
                eij = jw5(:,i) - jw5(:,j);
                jw = UU(:,ii - L + i) - UU(:,ii - L + j);%jw就是uij；
                tmp =1/(L^2*sigma2^2) * exp(-(eij^2)/2/sigma2^2) * eij * jw;
                w_hampel2 = w_hampel2 + mu6 * tmp;
            end
        end
    end
    
    %% robust5:Hampel 法

    en = zeros(L,1);
    for ii = L : LEN
        err_hampel3(mm,ii-L + 1) = (w_hampel3-wo)'*(w_hampel3-wo);
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_hampel3' * UU(:,ii - L + jj);
        end 
        u = 0.6745 * en / median(abs(en - median(en)));
        %median(A)是求矩阵A每列的中位数，若个数为偶数，则为最中间两个值的平均值；
        %en - median(en)：en为10*1的列向量，该向量每个元素都减去其中位数,加绝对值会使每个元素都变为非负数；
        %u的分母最终是一个数，u为en的倍数（10*1列向量）
        aa = 0.5;
        bb = 2;
        cc = 4;
        w5 = zeros(L);
        jw5 = zeros(1,L);
        for jj = 1 : L
            if abs(u(jj)) <= aa
                w5(jj) = 1;
                jw5(:,jj) = en(jj);
            elseif (aa < abs(u(jj))) && (abs(u(jj)) <=bb)
                w5(jj) = aa / abs(u(jj));
                jw5(:,jj) = en(jj) * w5(jj);
            elseif (bb < abs(u(jj))) && (abs(u(jj)) <=cc)
                w5(jj) = (aa / abs(u(jj))) * ((cc-abs(u(jj)))/(cc-bb));
                jw5(:,jj) = en(jj) * w5(jj);
            elseif abs(u(jj)) > cc
                w5(jj) = 0;
                jw5(:,jj) = 0;
            end
        end
        for i = 1 : L
            for j = 1 : L   
                %这个权值是乘在了en上，不加权重eij=en(i) - en(j);
                eij = jw5(:,i) - jw5(:,j);
                jw = UU(:,ii - L + i) - UU(:,ii - L + j);%jw就是uij；
                tmp =1/(L^2*sigma3^2) * exp(-(eij^2)/2/sigma3^2) * eij * jw;
                w_hampel3 = w_hampel3 + mu7 * tmp;
            end
        end
    end
    
    %% robust5:Hampel 法

    en = zeros(L,1);
    for ii = L : LEN
        err_hampel4(mm,ii-L + 1) = (w_hampel4-wo)'*(w_hampel4-wo);
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_hampel4' * UU(:,ii - L + jj);
        end 
        u = 0.6745 * en / median(abs(en - median(en)));
        %median(A)是求矩阵A每列的中位数，若个数为偶数，则为最中间两个值的平均值；
        %en - median(en)：en为10*1的列向量，该向量每个元素都减去其中位数,加绝对值会使每个元素都变为非负数；
        %u的分母最终是一个数，u为en的倍数（10*1列向量）
        aa = 0.5;
        bb = 2;
        cc = 4;
        w5 = zeros(L);
        jw5 = zeros(1,L);
        for jj = 1 : L
            if abs(u(jj)) <= aa
                w5(jj) = 1;
                jw5(:,jj) = en(jj);
            elseif (aa < abs(u(jj))) && (abs(u(jj)) <=bb)
                w5(jj) = aa / abs(u(jj));
                jw5(:,jj) = en(jj) * w5(jj);
            elseif (bb < abs(u(jj))) && (abs(u(jj)) <=cc)
                w5(jj) = (aa / abs(u(jj))) * ((cc-abs(u(jj)))/(cc-bb));
                jw5(:,jj) = en(jj) * w5(jj);
            elseif abs(u(jj)) > cc
                w5(jj) = 0;
                jw5(:,jj) = 0;
            end
        end
        for i = 1 : L
            for j = 1 : L   
                %这个权值是乘在了en上，不加权重eij=en(i) - en(j);
                eij = jw5(:,i) - jw5(:,j);
                jw = UU(:,ii - L + i) - UU(:,ii - L + j);%jw就是uij；
                tmp =1/(L^2*sigma4^2) * exp(-(eij^2)/2/sigma4^2) * eij * jw;
                w_hampel4 = w_hampel4 + mu8 * tmp;
            end
        end
    end
 
     disp(mm)
end
toc

figure(1),hold on;
box on;
wid = 2.5;
plot(10*log10(mean(err_mee1)),'-.','color',[0.7451 0.7451 0.7451],'LineWidth',wid)
plot(10*log10(mean(err_mee2)),'-.g','LineWidth',wid)
plot(10*log10(mean(err_mee3)),'-.b','LineWidth',wid)
plot(10*log10(mean(err_mee4)),'-.c','LineWidth',wid)
plot(10*log10(mean(err_hampel1)),'-.r','LineWidth',wid)
plot(10*log10(mean(err_hampel2)),'-.k','LineWidth',wid)
plot(10*log10(mean(err_hampel3)),'-.m','LineWidth',wid)
plot(10*log10(mean(err_hampel4)),'-.','color',[0.7451 0.451 0.451],'LineWidth',wid)
h=legend('MEE1 (\eta=1.1, L=15, {\bf\sigma=5})','MEE2 ({\eta=3}, L=15, {\bf\sigma=10})','MEE3 ({\eta=11}, L=15, {\bf\sigma=20})','MEE4 ({\eta=600}, L=15, {\bf\sigma=100})','MMEE-Hample1 ({\eta=2},  L=15, {\bf\sigma=5}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4)','MMEE-Hample2 ({\eta=9}, L=15, {\bf\sigma=10}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4)','MMEE-Hample3 ({\eta=35}, L=15, {\bf\sigma=20}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4)','MMEE-Hample4 ({\eta=1000}, L=15, {\bf\sigma=100}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4)');
set(h,'FontName','Times New Roman','FontSize',15,'FontWeight','normal');
xlabel('Iteration');
ylabel('Weight error power');

E_MEE1 = mean(10* log10(mean(err_mee1(end/2:end))));
E_MEE2 = mean(10* log10(mean(err_mee2(end/2:end))));
E_MEE3 = mean(10* log10(mean(err_mee3(end/2:end))));
E_MEE4 = mean(10* log10(mean(err_mee4(end/2:end))));
E_Mm1 = mean(10* log10(mean(err_hampel1(end/2:end))));
E_Mm2 = mean(10* log10(mean(err_hampel2(end/2:end))));
E_Mm3 = mean(10* log10(mean(err_hampel3(end/2:end))));
E_Mm4 = mean(10* log10(mean(err_hampel4(end/2:end))));



% figure(1),hold on;
% box on;
% plot(10*log10(mean(err_lms)))
% plot(10*log10(mean(err_lms1)))
% plot(10*log10(mean(err_lms2)))
% plot(10*log10(mean(err_lms3)))
% plot(10*log10(mean(err_lms5)))
% legend('MEE-normal','1Andrew’s ','2Biweight','3Cauchy','5Hampel');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2),hold on;
% box on;
% plot(10*log10(mean(err_lms)))
% plot(10*log10(mean(err_lms6)))
% plot(10*log10(mean(err_lms7)))
% plot(10*log10(mean(err_lms8)))
% plot(10*log10(mean(err_lms9)))
% plot(10*log10(mean(err_lms10)))
% legend('MEE-normal','6Huber','7Logistic','8Median','9Talworth','10Welsch');

