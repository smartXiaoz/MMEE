%mee -lms
clear all;
close all;
LEN = 2500;
m = 5;
L=10;
sigma1 = 2;
sigma2 = 10;
sigma3 = 20;
sigma4 = 2.5;
mu0 = 0.045;
mu1 = 0.6;
mu2 = 3;
mu3 = 11;
mu4 = 400;
mu5 = 2;
mu6 = 9;
mu7 = 35;
mu8 = 0.4;
tic
for mm = 1 : 100
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
    w_lms0 = w_lms;
   
    
  
    
     
    
    %% robust5:Hampel 法
    Dm2 = zeros(m,m);
    Dm3 = zeros(m,1);
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
            Dm2 = Dm2 + ul*ww*(P-Q)*ww*ul';
            Dm3 = Dm3 + ul*ww*(P-Q)*ww*vvll;
        end
    end
    
    LEN1 = LEN - 500;

    
    
    eep = ((Dm3'*Dm3)/LEN1);

    yita =  1/(L^2*sigma4^2);
    mu11 = yita*mu8;
    aaa = mu11^2*eep;
    eeu = mean(diag(Dm2/LEN1));
    eeccc = 1-(1-mu11*eeu)^2;
    stdeeem(mm) = aaa/eeccc;
    
     disp(mm)
    
    
end
toc






figure(1),hold on;
box on;
wid = 1.5;

% plot(10*log10(mean(err_mee4)),'-.c','LineWidth',wid)
% % plot(10*log10(mean(err_hampel1)),'-.r','LineWidth',wid)
% % plot(10*log10(mean(err_hampel2)),'-.k','LineWidth',wid)
% % plot(10*log10(mean(err_hampel3)),'-.m','LineWidth',wid)
plot(10*log10(mean(err_hampel4)),'-.r','LineWidth',wid)
plot(10*log10(ones(1,LEN)*mean(stdeeem)),'-.r','LineWidth',wid)
