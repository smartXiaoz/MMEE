clear all; close all;
%% MEE-AF-V10
L = 2000;
p = 5;
q = 8;
lamda1 = 0.997;%RMEE 
lamda2 = 0.996;%RMEE
lamda3 = 0.990;%RMEE
lamda4 = 0.985;%RMEE
lamda5 = 0.970;%RMEE

% lamda2 = 1;%RMEE for speech signal

sigma1= 0.2;%%RMEE
sigma2 = 0.5;%%RMEE
sigma3 = 1;%%RMEE
sigma4 = 2;%%RMEE
sigma5 = 3;%%RMEE
tic
for mm = 1:2
    vv = randn(1,L) * 0.1;
%     vv = rand(1,L) -0.5;
    v1=randn(1,L)*0.1; v2=randn(1,L)*5;
    rp=rand(1,L);
    %    vv = v1 + (rp>0.95).*v2;
    vv = (rp<=0.95).*v1 + (rp>0.95).*v2;
    vG = exp(-vv.^2/2/sigma1^2).*vv;
    %      vv = exp(-vv.^2/2/sigma1^2).*vv;
    
   wo1 = randn(p,1);wo2 = randn(p,1);wo3 = randn(p,1);
%       wo1 = [0,0,0.9,0,0,0,0.2,0,0,0]';
%     wo = [ kron(wo1, ones(1,L/3)) kron(wo2, ones(1,L/3)) kron(wo3, ones(1,L/3))];
     wo = [ kron(wo1, ones(1,L)) ];
  %   wo = [ kron(randn(p,1), ones(1,L/2)),  kron(randn(p,1), ones(1,L/2)) ];
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
    w_M_RLS1 =w_LMS;%%MEE
    w_M_RLS2 =w_LMS;%%MEE
    w_M_RLS3 =w_LMS;%%MEE
    w_M_RLS4 =w_LMS;%%MEE
    w_M_RLS5 =w_LMS;%%MEE
    w_M_RLS6 =w_LMS;%%MEE
    
    %% RLS MEE, algorithm 1
    PL = eye(p)*1; 
    eL = zeros(q, 1); 
    PHI = zeros(q, q);
    for ii = q : L
        Err_RLS_MEE1(mm,ii) = (wo(:,ii)  - w_M_RLS1)' * (wo(:,ii)  - w_M_RLS1);
        for jj = 1 : q
            eL(jj) = dd(ii - q + jj) - w_M_RLS1' * uu(: , ii - q + jj);
        end
        u0 = uu(:,ii);
        e0 = eL(q);
        phi0 = 0;
        for kk = 2:q
            ek = eL(q-kk+1);
            phi0 = phi0 + lamda1^(kk-1) * exp(- (e0-ek).^2/sigma1^2/2);
        end
        KL = PL * u0 * inv ( lamda1^2 / phi0 + u0' * PL * u0 );
        PL = 1/lamda1^2 * (PL - KL * u0' * PL);
        w_M_RLS1 = w_M_RLS1 + KL * e0;        
    end
     %% RLS MEE, algorithm 2
    PL = eye(p)*1; 
    eL = zeros(q, 1); 
    PHI = zeros(q, q);
    for ii = q : L
        Err_RLS_MEE2(mm,ii) = (wo(:,ii)  - w_M_RLS2)' * (wo(:,ii)  - w_M_RLS2);
        for jj = 1 : q
            eL(jj) = dd(ii - q + jj) - w_M_RLS2' * uu(: , ii - q + jj);
        end
        u0 = uu(:,ii);
        e0 = eL(q);
        phi0 = 0;
        for kk = 2:q
            ek = eL(q-kk+1);
            phi0 = phi0 + lamda2^(kk-1) * exp(- (e0-ek).^2/sigma2^2/2);
        end
        KL = PL * u0 * inv ( lamda2^2 / phi0 + u0' * PL * u0 );
        PL = 1/lamda2^2 * (PL - KL * u0' * PL);
        w_M_RLS2 = w_M_RLS2 + KL * e0;        
    end
     %% RLS MEE, algorithm 3
    PL = eye(p)*1; 
    eL = zeros(q, 1); 
    PHI = zeros(q, q);
    for ii = q : L
        Err_RLS_MEE3(mm,ii) = (wo(:,ii)  - w_M_RLS3)' * (wo(:,ii)  - w_M_RLS3);
        for jj = 1 : q
            eL(jj) = dd(ii - q + jj) - w_M_RLS3' * uu(: , ii - q + jj);
        end
        u0 = uu(:,ii);
        e0 = eL(q);
        phi0 = 0;
        for kk = 2:q
            ek = eL(q-kk+1);
            phi0 = phi0 + lamda3^(kk-1) * exp(- (e0-ek).^2/sigma3^2/2);
        end
        KL = PL * u0 * inv ( lamda3^2 / phi0 + u0' * PL * u0 );
        PL = 1/lamda3^2 * (PL - KL * u0' * PL);
        w_M_RLS3 = w_M_RLS3 + KL * e0;        
    end
    
    %% RLS MEE, algorithm 4
    PL = eye(p)*1; 
    eL = zeros(q, 1); 
    PHI = zeros(q, q);
    for ii = q : L
        Err_RLS_MEE4(mm,ii) = (wo(:,ii)  - w_M_RLS4)' * (wo(:,ii)  - w_M_RLS4);
        for jj = 1 : q
            eL(jj) = dd(ii - q + jj) - w_M_RLS4' * uu(: , ii - q + jj);
        end
        u0 = uu(:,ii);
        e0 = eL(q);
        phi0 = 0;
        for kk = 2:q
            ek = eL(q-kk+1);
            phi0 = phi0 + lamda4^(kk-1) * exp(- (e0-ek).^2/sigma4^2/2);
        end
        KL = PL * u0 * inv ( lamda4^2 / phi0 + u0' * PL * u0 );
        PL = 1/lamda4^2 * (PL - KL * u0' * PL);
        w_M_RLS4 = w_M_RLS4 + KL * e0;        
    end
     %% RLS MEE, algorithm 5
    PL = eye(p)*1; 
    eL = zeros(q, 1); 
    PHI = zeros(q, q);
    for ii = q : L
        Err_RLS_MEE5(mm,ii) = (wo(:,ii)  - w_M_RLS5)' * (wo(:,ii)  - w_M_RLS5);
        for jj = 1 : q
            eL(jj) = dd(ii - q + jj) - w_M_RLS5' * uu(: , ii - q + jj);
        end
        u0 = uu(:,ii);
        e0 = eL(q);
        phi0 = 0;
        for kk = 2:q
            ek = eL(q-kk+1);
            phi0 = phi0 + lamda5^(kk-1) * exp(- (e0-ek).^2/sigma5^2/2);
        end
        KL = PL * u0 * inv ( lamda5^2 / phi0 + u0' * PL * u0 );
        PL = 1/lamda5^2 * (PL - KL * u0' * PL);
        w_M_RLS5 = w_M_RLS5 + KL * e0;        
    end
    
     %% RLS MEE, algorithm 6
    PL = eye(p)*1; 
    eL = zeros(q, 1); 
    PHI = zeros(q, q);
    for ii = q : L
        
        Err_RLS_MEE6(mm,ii) = (wo(:,ii)  - w_M_RLS6)' * (wo(:,ii)  - w_M_RLS6);
        for jj = 1 : q
            eL(jj) = dd(ii - q + jj) - w_M_RLS6' * uu(: , ii - q + jj);
        end
        u0 = uu(:,ii);
        e0 = eL(q);
        phi0 = 0;
        en = eL;
        
        u = 0.6745 * en / median(abs(en - median(en)));
        %median(A)是求矩阵A每列的中位数，若个数为偶数，则为最中间两个值的平均值；
        %en - median(en)：en为10*1的列向量，该向量每个元素都减去其中位数,加绝对值会使每个元素都变为非负数；
        %u的分母最终是一个数，u为en的倍数（10*1列向量）
        aa = 0.2;
        bb = 1;
        cc = 1.5;
        w5 = zeros(q);
        jw5 = zeros(1,q);
        for jj = 1 : q
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
        
        
        for kk = 2:q
            if ii > 600 % 开始不使用鲁棒
            ek = jw5(:,q-kk+1);
            e0 = jw5(:,q);
            else
            ek = eL(q-kk+1);
            e0 = eL(q);
            end
            phi0 = phi0 + lamda1^(kk-1) * exp(- (e0-ek).^2/sigma1^2/2);
        end
        KL = PL * u0 * inv ( lamda1^2 / phi0 + u0' * PL * u0 );
        PL = 1/lamda1^2 * (PL - KL * u0' * PL);
        w_M_RLS6 = w_M_RLS6 + KL * e0;        
    end
    
    %% RLS-MEE, steady-state mean square performance
    Gsig =0; %phi0
    phi01 = 0;  phi02 = 0;  phi03 = 0; phi04 = 0; phi05 = 0;
    Gsig1 =0; Gsig2 =0; Gsig3 =0; Gsig4 =0; Gsig5 =0;
    v0 = vv(q: end);
    for ii = 2:q
        phi01 = phi01 + lamda1^(ii-1);
        phi02 = phi02 + lamda2^(ii-1);
        phi03 = phi03 + lamda3^(ii-1);
        phi04 = phi04 + lamda4^(ii-1);
        phi05 = phi05 + lamda5^(ii-1);
    end
    for ii = 2:q
        vi = vv(q-ii+1 : end-ii+1);
        Gsig1 = Gsig1 + lamda1^(ii-1) * exp(- (v0-vi).^2/sigma1^2/2);
        Gsig2 = Gsig2 + lamda2^(ii-1) * exp(- (v0-vi).^2/sigma2^2/2);
        Gsig3 = Gsig3 + lamda3^(ii-1) * exp(- (v0-vi).^2/sigma3^2/2);
        Gsig4 = Gsig4 + lamda4^(ii-1) * exp(- (v0-vi).^2/sigma4^2/2);
        Gsig5 = Gsig5 + lamda5^(ii-1) * exp(- (v0-vi).^2/sigma5^2/2);
    end
    tmp1 = (v0.*Gsig1) ;
    tmp2 = (v0.*Gsig2) ;
    tmp3 = (v0.*Gsig3) ;
    tmp4 = (v0.*Gsig4) ;
    tmp5 = (v0.*Gsig5) ;
    L1 = length(tmp1);
     
    a =(1-lamda1^2)/(1+lamda1^2) * p / mean(diag(uu*uu'/L)) ...
        * tmp1*tmp1'/L1 / mean(Gsig1)^2; %phi01^2; %
    Err_TH_MEE1(mm) = a;
    
    a =(1-lamda2^2)/(1+lamda2^2) * p / mean(diag(uu*uu'/L)) ...
        * tmp2*tmp2'/L1 / mean(Gsig2)^2; %phi02^2; %
    Err_TH_MEE2(mm) = a;
    
    a =(1-lamda3^2)/(1+lamda3^2) * p / mean(diag(uu*uu'/L)) ...
        * tmp3*tmp3'/L1 / mean(Gsig3)^2; %phi03^2; %
    Err_TH_MEE3(mm) = a;
    
    a =(1-lamda4^2)/(1+lamda4^2) * p / mean(diag(uu*uu'/L)) ...
        * tmp4*tmp4'/L1 / mean(Gsig4)^2; %phi04^2; %
    Err_TH_MEE4(mm) = a;
    
    a =(1-lamda5^2)/(1+lamda5^2) * p / mean(diag(uu*uu'/L)) ...
        * tmp5*tmp5'/L1 / mean(Gsig5)^2; %phi05^2; %
    Err_TH_MEE5(mm) = a;
    
end



toc
figure,hold on;
plot(10* log10(mean(Err_RLS_MEE1)),'r'),
plot(10* log10(mean(Err_RLS_MEE2)),'g'),
plot(10* log10(mean(Err_RLS_MEE3)),'b');
plot(10* log10(mean(Err_RLS_MEE4)),'c'),
plot(10* log10(mean(Err_RLS_MEE6)),'k');
plot(10*log10(ones(1,L)*mean(Err_TH_MEE1)),'r--')
plot(10*log10(ones(1,L)*mean(Err_TH_MEE2)),'g--')
plot(10*log10(ones(1,L)*mean(Err_TH_MEE3)),'b--')
plot(10*log10(ones(1,L)*mean(Err_TH_MEE4)),'c--')
plot(10*log10(ones(1,L)*mean(Err_TH_MEE5)),'k--')
legend('RMEE1','RMEE2','RMEE3','RMEE4','RMEE5', ...
'TH-RMEE1', 'TH-RMEE2', 'TH-RMEE3', 'TH-RMEE4', 'TH-RMEE5');
xlabel('Iterations');ylabel('MSD');
% legend('LMS','RLS','RMC','RMEE');
[10*log10(mean(Err_TH_MEE1))
10*log10(mean(Err_TH_MEE2))
10*log10(mean(Err_TH_MEE3))
10*log10(mean(Err_TH_MEE4))
10*log10(mean(Err_TH_MEE5))]'