%mee -lms
clear all;
close all;
LEN = 2000;
m = 5;
L=10;
% sigma1 = 2;
% sigma2 = 10;
% sigma3 = 20;
sigma4 = 2.5;
% mu0 = 0.045;
% mu1 = 0.2;
% mu2 = 3;
% mu3 = 11;
% mu4 = 400;
% mu5 = 2;
% mu6 = 9;
% mu7 = 35;
% mu8 = 0.1;
tic
t_mean = 0;
f_mean = 0;
mmii = 0;
etamm = [0.1];
% etamm = [0.06 0.07 0.08 0.09 0.1 0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 0.20 0.21 0.22];
% etamm = [0.06 0.16 0.26 0.36 0.46 0.56 0.66 0.76 0.86];
for mmu = etamm
    mmii = mmii + 1;
    
for mm = 1 : 2
    mu8 = mmu;
    wo = randn(m, 1);%初始w0为30*1
    w_hampel4 = randn(m, 1);
    %vv = randn(1,LEN) * 0.1;
    e_greedy = rand(1, LEN);%均匀分布
    vv=randn(1,LEN)*0.1;
    v1 = randn(1,LEN) * 0.1;
    v2 = randn(1,LEN) * 10;
    VV = (e_greedy <= 0.99).*v1 + (e_greedy >0.99).*v2; %non-gaussian-noise
%     VV = v1;
   
    UU = randn(m, LEN);%uu为输入，理论上为30*1，但取3000组数据；
    alpha_noise = alpha_stable_noise(1.1,0.1,0,0,LEN);%alpha=1.1 gamma=0.1
    varNoise = 1;
    alpha_noise = sqrt(varNoise) * alpha_noise;
    for ii = 1 : LEN
        DD(ii) = wo' * UU(:,ii) + VV(ii);%产生1*3000的行向量，也即是输出
%         DD(ii) = wo' * UU(:,ii) +alpha_noise(ii);%产生1*3000的行向量，也即是输出
    end 
%% robust5:Hampel 法
    Dm2 = zeros(m,m);
    
    Dm3 = zeros(m,1);
    Dm4 = 0;
    Dm5 = 0;
    DM6 = 0;
    DW0 = 0;
    T_1 = 0;
    T_2 = 0;
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
                uul(:,jj) = ul(:,jj)*w5(jj);
%                 vvll(jj) = vvll(jj)*w5(jj)*w5(jj);
            elseif (bb < abs(u(jj))) && (abs(u(jj)) <=cc)
                w5(jj) = (aa / abs(u(jj))) * ((cc-abs(u(jj)))/(cc-bb));
                jw5(:,jj) = en(jj) * w5(jj);
                uul(:,jj) = ul(:,jj)*w5(jj);
%                 vvll(jj) = vvll(jj)*w5(jj)*w5(jj);
            elseif abs(u(jj)) > cc
                w5(jj) = 0;
                jw5(:,jj) = 0;
                uul(:,jj) = ul(:,jj)*w5(jj);
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
%                 jw = uul(:,i) - uul(:,j);
                g = exp(-(eij^2 )/ 2 / sigma4^2);
                Q(i, j) = g;
                P(i, i) = P(i, i) + g;
                tmp =1/(L^2*sigma4^2) * g * eij * jw;
                w_hampel4 = w_hampel4 + mu8 * tmp;
            end
        end
        if ii > 1000
            ww = diag(w5);
            DW0 = DW0 + ww;
            Dm2 = Dm2 + ul*ww*(P-Q)*ww*ul';
            Dm3 = Dm3 + ul*ww*(P-Q)*ww*vvll;
            Dm4 = Dm4 + (ul*ww*(P-Q)*ww*vvll)'*(ul*ww*(P-Q)*ww*vvll);
            Dm5 = Dm5 + ul'*ul;
            DM6 = DM6 + vvll'*ww*(P-Q)*ww*ww*(P-Q)*ww*vvll;
            T_1 = T_1 + Q;
            T_2 = T_2 + P-Q;
        end
    end
    D6 = 0;
    Mvv = mhampel(VV,LEN);
    for i = 1 : L-1
        for j = 1000  : LEN 
           vij = Mvv(j) - Mvv(j-i);
           D6 = D6 + exp(-(vij^2 )/ 2 / sigma4^2);
        end
    end
    LEN1 = LEN-1000;
%     for i = 2 : LEN       
%            vij = VV(i) - VV(i-1);
%            D6 = D6 + exp(-(vij^2 )/ 2 / sigma1^2);       
%     end
    D6 = D6/LEN1;%这里算的phi的均值,UL'UL = m  
    
%     Dm4 = mean(diag(Dm4/LEN1));
%     D6 = Dm4;
%     eep = ((Dm3'*Dm3)/LEN1);
    
%     eep = Dm4*L/LEN1;
    tmp1 = mean(diag(DW0/LEN1));
    eep = DM6/LEN1*m*L*tmp1^2;
%     eep = mean(diag(DM6/LEN1))*mean(diag(Dm5/LEN1))*m*tmp1^2;
%     eep = mean(diag(Dm4/LEN1));
%     eep = (D6*tmp1^2)^2 * var(VV) * m;
%     eep = (D6*tmp1^8)^2 * var(VV) * m;
    yita =  1/(L^2*sigma4^2);
    mu11 = yita*mu8;
    aaa = mu11^2*eep;       
    eeu = D6*tmp1^2*L;    
    eeccc = 1-(1-mu11*2*eeu)^2;
    stdeeem(mm) = 4*aaa/eeccc;
%     eep = ((Dm3'*Dm3)/LEN1);
%     yita =  1/(L^2*sigma4^2);
%     mu11 = yita*mu8;
%     aaa = mu11^2*eep;
%     eeu = mean(diag(Dm2/LEN1));
%     eeccc = 1-(1-mu11*eeu)^2;
%     stdeeem(mm) = aaa/eeccc;
    f_mean(mm,mmii) = mean(mean(err_hampel4(:,end-1000:end)));
    t_mean(mm,mmii) = mean(stdeeem);
    disp(mm)
end  
mean(f_mean(:,mmii))
mean(t_mean(:,mmii))
end  
figure,hold on;
box on;
wid = 1.5;
plot(10*log10(mean(err_hampel4)),'-.r','LineWidth',wid)
plot(10*log10(ones(1,LEN)*mean(stdeeem)),'-.r','LineWidth',wid)

figure,hold on;
box on;
wid = 2.5;
% plot(etamm,10*log10(f_mean),'-.r','LineWidth',wid)
% plot(etamm,10*log10(t_mean),'-.b','LineWidth',wid)
plot(etamm,mean(f_mean),'-r','LineWidth',wid)
plot(etamm,mean(t_mean),':bO','LineWidth',wid)
legend('simulation','theroy');
xlabel('\eta');
ylabel('MSD');
%     