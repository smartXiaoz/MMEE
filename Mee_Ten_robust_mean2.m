%mee -lms
clear all;
close all;
LEN = 10000;
m = 5;
L=10;
sigma1 = 2;
sigma2 = 10;
sigma3 = 20;
sigma4 = 2;
mu0 = 0.045;
mu1 = 0.2;
mu2 = 3;
mu3 = 11;
mu4 = 400;
mu5 = 2;
mu6 = 9;
mu7 = 35;
mu8 = 0.1;
tic
for mm = 1 : 20
    e_greedy = rand(1, LEN);%均匀分布
    %vv = randn(1,LEN) * 0.1;
    vv=randn(1,LEN)*0.1;
    v1 = randn(1,LEN) * 0.1;
    v2 = randn(1,LEN) * 10;
    VV = (e_greedy <= 1.9).*v1 + (e_greedy >1.9).*v2; %non-gaussian-noise
%     VV = v1;
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
    w_lms0 = w_lms;
   
    
   
    sum = 0;
    sum2 = 0;
    D1 = zeros(m,m);
    D2 = zeros(m,m);
    D3 = zeros(m,1);
    D4 = 0;
    D5 = 0;
    D7 = 0;
    D8 = 0;
    D9 = 0;
    for ii = L : LEN
         %%MEE
        en = zeros(L,1);%L=10；   
        ul = zeros(m,L);
        vvll = zeros(L,1);
        err_mee1(mm,ii-L + 1) = (w_mee1-wo)'*(w_mee1-wo);%mm为最外层循环，从1-200；最终就要画出此图像；
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_mee1' * UU(:,ii - L + jj);
            ul(:,jj) = UU(:,ii - L + jj);
            %en这个列向量始终存储10个误差1-10，2-11，3-12......
            vvll(jj) = VV(:,ii - L + jj);
            
        end 
        
        sum = sum + ul*ul';
        sum2 = sum2 + ul'*ul;
        P = zeros(L, L);
        Pv = zeros(L, L);
        Q = zeros(L, L);
        for i = 1 : L
            for j = 1 : L
                uij = UU(:,ii - L + i) - UU(:,ii - L + j);
                eij = en(i) - en(j);
                vij = vvll(i) - vvll(j);
                g = exp(-(eij^2 )/ 2 / sigma1^2);
                g1 = exp(-(vij^2 )/ 2 / sigma1^2);
                Q(i, j) = g;
                P(i, i) = P(i, i) + g;
                Pv(i, i) = Pv(i, i) + g1;
                tmp =  1/(L^2*sigma1^2)*g * eij * uij;
                w_mee1 = w_mee1 + mu1 * tmp;%mul=0.2;sigma=2;
            end           
        end
        if ii > 0
            D1 = D1 + ul*(P-eye(L))*ul';
            D2 = D2 + ul*(P-Q)*ul';
            D3 = D3 + ul*(P-Q)*vvll;
%             D4 = D4 + (P-eye(L));
%             D5 = D5 + ul*(Pv-eye(L))*vvll;
%             D7 = D7 + ul'*ul;
            D4 = D4 + Q;
            D5 = D5 + (P-Q)*(P-Q);
            D7 = D7 + vvll'*(P-Q)*(P-Q)*vvll;
            D8 = D8 + (eye(m)-mu1*ul*(P-Q)*ul')*(eye(m)-mu1*ul*(P-Q)*ul');
            D9 = D9 + vvll'*(Pv-Q)*ul'*ul*(Pv-Q)*vvll;
            
%             D5 = D5 + P;
        end
        
    end
    
     
    
    %% robust5:Hampel 法
    Dm2 = zeros(m,m);
    
    Dm3 = zeros(m,1);
    Dm4 = 0;
    DW0 = 0;
    en = zeros(L,1);
    for ii = L : LEN
        en = zeros(L,1);%L=10；   
        ul = zeros(m,L);
        uul = zeros(m,L);
        vvll = zeros(L,1);
        ww = zeros(L,L);
        err_hampel4(mm,ii-L + 1) = (w_hampel4-wo)'*(w_hampel4-wo);
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_hampel4' * UU(:,ii - L + jj);
            ul(:,jj) = UU(:,ii - L + jj);
            vvll(jj) = VV(:,ii - L + jj);
        end 
        u = 0.6745 * en / median(abs(en - median(en)));
        %median(A)是求矩阵A每列的中位数，若个数为偶数，则为最中间两个值的平均值；
        %en - median(en)：en为10*1的列向量，该向量每个元素都减去其中位数,加绝对值会使每个元素都变为非负数；
        %u的分母最终是一个数，u为en的倍数（10*1列向量）
        aa = 2;
        bb = 4;
        cc = 8;
        w5 = zeros(L,1);
        jw5 = zeros(1,L);
        for jj = 1 : L
            if abs(u(jj)) <= aa
                w5(jj) = 1;
                jw5(:,jj) = en(jj);
                uul(:,jj) = ul(:,jj);
                
            elseif (aa < abs(u(jj))) && (abs(u(jj)) <=bb)
                w5(jj) = aa / abs(u(jj));
                jw5(:,jj) = en(jj) * w5(jj);
                uul(:,jj) = ul(:,jj)*w5(jj)*w5(jj);
%                 vvll(jj) = vvll(jj)*w5(jj)*w5(jj);
            elseif (bb < abs(u(jj))) && (abs(u(jj)) <=cc)
                w5(jj) = (aa / abs(u(jj))) * ((cc-abs(u(jj)))/(cc-bb));
                jw5(:,jj) = en(jj) * w5(jj);
                uul(:,jj) = ul(:,jj)*w5(jj)*w5(jj);
%                 vvll(jj) = vvll(jj)*w5(jj)*w5(jj);
            elseif abs(u(jj)) > cc
                w5(jj) = 0;
                jw5(:,jj) = 0;
                uul(:,jj) = ul(:,jj)*w5(jj)*w5(jj);
%                 vvll(jj) = vvll(jj)*w5(jj)*w5(jj);
            end
        end
        P = zeros(L, L);
        Q = zeros(L, L);
        for i = 1 : L
            for j = 1 : L   
                %这个权值是乘在了en上，不加权重eij=en(i) - en(j);
                eij = jw5(:,i) - jw5(:,j);
                jw = UU(:,ii - L + i) - UU(:,ii - L + j);%jw就是uij；
                
                g = exp(-(eij^2 )/ 2 / sigma4^2);
                Q(i, j) = g;
                P(i, i) = P(i, i) + g;
                tmp =1/(L^2*sigma4^2) * g * eij * jw;
                w_hampel4 = w_hampel4 + mu8 * tmp;
            end
        end
        if ii > 500
            ww = diag(w5);
            DW0 = DW0 + ww;
            Dm2 = Dm2 + ul*ww*(P-Q)*ww*ul';
            Dm3 = Dm3 + ul*ww*(P-Q)*ww*vvll;
            Dm4 = Dm4 + ww*(P-Q)*ww*vvll;
        end
    end
    
    D6 = 0;
    for i = 1 : L-1
        for j = L : LEN 
           vij = VV(j) - VV(j-i);
           D6 = D6 + exp(-(vij^2 )/ 2 / sigma1^2);
        end
    end
%     for i = 2 : LEN       
%            vij = VV(i) - VV(i-1);
%            D6 = D6 + exp(-(vij^2 )/ 2 / sigma1^2);       
%     end
    D6 = D6/LEN;%这里算的phi的均值,UL'UL = m  
    
    D10 = 0;
    for i = 1 : L-1
        tmp = 0;
        tmp1 = 0;
        for j = L : LEN 
           vij = VV(j) - VV(j-i);
           tmp = tmp + exp(-(vij^2 )/ 2 / sigma1^2);
           tmp1 = tmp1 + (exp(-(vij^2 )/ 2 / sigma1^2))^2;
        end
        D10 = D10 + tmp^2 + tmp1;
    end
    D10 = D10/LEN;%这里算的phi的均值,UL'UL = m  
    
%     DM6 = 0;
%     u = 0.6745 * VV / median(abs(en - median(VV)));
%         %median(A)是求矩阵A每列的中位数，若个数为偶数，则为最中间两个值的平均值；
%         %en - median(en)：en为10*1的列向量，该向量每个元素都减去其中位数,加绝对值会使每个元素都变为非负数；
%         %u的分母最终是一个数，u为en的倍数（10*1列向量）
%         aa = 2;
%         bb = 4;
%         cc = 8;
%         w5 = zeros(L,1);
%         jw5 = zeros(1,L);
%         for jj = 1 : L
%             if abs(u(jj)) <= aa
%                 w5(jj) = 1;
%                 jw5(:,jj) = en(jj);
%                 uul(:,jj) = ul(:,jj);
%                 
%             elseif (aa < abs(u(jj))) && (abs(u(jj)) <=bb)
%                 w5(jj) = aa / abs(u(jj));
%                 jw5(:,jj) = en(jj) * w5(jj);
%                 uul(:,jj) = ul(:,jj)*w5(jj)*w5(jj);
% %                 vvll(jj) = vvll(jj)*w5(jj)*w5(jj);
%             elseif (bb < abs(u(jj))) && (abs(u(jj)) <=cc)
%                 w5(jj) = (aa / abs(u(jj))) * ((cc-abs(u(jj)))/(cc-bb));
%                 jw5(:,jj) = en(jj) * w5(jj);
%                 uul(:,jj) = ul(:,jj)*w5(jj)*w5(jj);
% %                 vvll(jj) = vvll(jj)*w5(jj)*w5(jj);
%             elseif abs(u(jj)) > cc
%                 w5(jj) = 0;
%                 jw5(:,jj) = 0;
%                 uul(:,jj) = ul(:,jj)*w5(jj)*w5(jj);
% %                 vvll(jj) = vvll(jj)*w5(jj)*w5(jj);
%             end
%         end
%     for i = 1 : L-1
%         for j = L : LEN 
%            vij = VV(j) - VV(j-i);
%            DM6 = DM6 + exp(-(vij^2 )/ 2 / sigma1^2);
%         end
%     end
%     DM6 = DM6/LEN;%这里算的phi的均值,UL'UL = m 
    
    
    LEN1 = LEN - 500;
    eep = ((D3'*D3)/LEN1);
    
%     eep = D9/LEN1*m;
    eep = D6^2 * var(VV) * m;
    yita =  1/(L^2*sigma1^2);
    mu11 = yita*mu1;
    aaa = mu11^2*eep;
    eeu = mean(diag(D2/LEN1));
    eeu = D6*L;
    eeccc = 1-(1-mu11*eeu)^2;
%     eeccc = 1-(mean(diag(D2/LEN1)));
    stdeee(mm) = aaa/eeccc;
       
    
    
    
    eep = ((Dm3'*Dm3)/LEN1);
%     eep = (Dm4'*Dm4*m)/LEN1;
    tmp1 = mean(diag(DW0/LEN1));
    eep = (D6*tmp1^8)^2 * var(VV) * m;
    yita =  1/(L^2*sigma4^2);
    mu11 = yita*mu8;
    aaa = mu11^2*eep;
    
    
%     eeu = D6*L*tmp1^2;
    
    eeccc = 1-(1-mu11*eeu)^2;
    stdeeem(mm) = aaa/eeccc;
    
    disp(mm)
    
    
end
toc





mean(stdeee)
mean(mean(err_mee1))

figure(1),hold on;
box on;
wid = 1.5;
% 
% plot(10*log10(mean(err_lms0)),'-.b','LineWidth',wid)
% % plot(10*log10(mean(err_mee1)),'-.','color',[0.7451 0.7451 0.7451],'LineWidth',wid)
plot(10*log10(mean(err_mee1)),'-.b','LineWidth',wid)
plot(10*log10(ones(1,LEN)*mean(stdeee)),'-.b','LineWidth',wid)
% plot(10*log10(mean(err_mee4)),'-.c','LineWidth',wid)
% % plot(10*log10(mean(err_hampel1)),'-.r','LineWidth',wid)
% % plot(10*log10(mean(err_hampel2)),'-.k','LineWidth',wid)
% % plot(10*log10(mean(err_hampel3)),'-.m','LineWidth',wid)
plot(10*log10(mean(err_hampel4)),'-.r','LineWidth',wid)
plot(10*log10(ones(1,LEN)*mean(stdeeem)),'-.r','LineWidth',wid)
% plot(10*log10(ones(1,LEN)*mean(stdeee)),'-.b','LineWidth',wid)
% %h=legend('LMS','MEE1 (\eta=1.1, L=15, {\bf\sigma=5})','MEE2 ({\eta=3}, L=15, {\bf\sigma=10})','MEE3 ({\eta=11}, L=15, {\bf\sigma=20})','MEE4 ({\eta=600}, L=15, {\bf\sigma=100})','MMEE-Hample1 ({\eta=2},  L=15, {\bf\sigma=5}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4)','MMEE-Hample2 ({\eta=9}, L=15, {\bf\sigma=10}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4)','MMEE-Hample3 ({\eta=35}, L=15, {\bf\sigma=20}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4)','MMEE-Hample4 ({\eta=1000}, L=15, {\bf\sigma=100}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4)');
% h=legend('LMS (\eta=0.045)','MEE (\eta=0.7, L=15, {\sigma=2})','MEE ({\eta=400}, L=15, {\sigma=100})','MMEE-Hample ({\eta=400}, L=15, {\sigma=100}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4)');
% set(h,'FontName','Times New Roman','FontSize',15,'FontWeight','normal');
% xlabel('Iteration');
% ylabel('Weight error power');

% E_MEE1 = mean(10* log10(mean(err_mee1(end/2:end))));
% E_MEE2 = mean(10* log10(mean(err_mee2(end/2:end))));
% E_MEE3 = mean(10* log10(mean(err_mee3(end/2:end))));
% E_MEE4 = mean(10* log10(mean(err_mee4(end/2:end))));
% E_Mm1 = mean(10* log10(mean(err_hampel1(end/2:end))));
% E_Mm2 = mean(10* log10(mean(err_hampel2(end/2:end))));
% E_Mm3 = mean(10* log10(mean(err_hampel3(end/2:end))));
% E_Mm4 = mean(10* log10(mean(err_hampel4(end/2:end))));



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
