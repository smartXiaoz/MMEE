%mee -lms
clear all;
close all;
LEN = 1500;
m = 5;
L = 15;

mu_lms = 0.02;
mu_lmm = 0.02;
mu_gmcc = 0.019;
mu_mee = 0.5;
mu_andrews = 0.3;
mu_huber = 0.4;
mu_hampel = 0.3;

%case3
% mu_lms = 0.02;
% mu_lmm = 0.02;
% mu_gmcc = 0.019;
% mu_mee = 0.4;
% mu_andrews = 0.3;
% mu_huber = 0.3;
% mu_hampel = 0.3;

sigma = 2; 
tic
Epoch = 100;
err_lms0 = zeros(Epoch, LEN);
err_lms11 = zeros(Epoch, LEN);
Err_GMCC_LMS = zeros(Epoch, LEN);
err_lms = zeros(Epoch, LEN);
err_lms1 = zeros(Epoch, LEN);
err_lms6 = zeros(Epoch, LEN);
err_lms5 = zeros(Epoch, LEN);
for mm = 1 :Epoch
% err_lms0 = zeros(Epoch, LEN);
% err_lms11 = zeros(Epoch, LEN);
% Err_GMCC_LMS = zeros(Epoch, LEN);
% err_lms = zeros(Epoch, LEN);
% err_lms1 = zeros(Epoch, LEN);
% err_lms6 = zeros(Epoch, LEN);
% err_lms5 = zeros(Epoch, LEN);
    e_greedy = rand(1, LEN);%���ȷֲ�
    %vv = randn(1,LEN) * 0.1;
    wo = randn(m, 1);%��ʼw0Ϊ30*1
    UU = randn(m, LEN);%uuΪ���룬������Ϊ30*1����ȡ3000�����ݣ�
    DD = zeros(LEN,1);
    %% case1
    vv=randn(1,LEN)*0.1;
    v1 = randn(1,LEN) * 0.1;
    v2 = randn(1,LEN) * 10;
    VV = (e_greedy <= 0.90).*v1 + (e_greedy >0.90).*v2; %non-gaussian-noise
    
    %% case2

%     alpha_noise = alpha_stable_noise(1.2,0.1,0,0,LEN);%alpha=1.1 gamma=0.1
%     varNoise = 1;
%     VV = sqrt(varNoise) * alpha_noise;

    %% case3 chi-r
%     VV = chi2rnd(1,[1, LEN]) * 0.1;


    %% case4 t_dis
    VV = trnd(1.5,[1, LEN]) * 0.1;
    for ii = 1 : LEN
        DD(ii) = wo' * UU(:,ii) + VV(ii);%����1*3000����������Ҳ�������
%         DD(ii) = wo' * UU(:,ii) +alpha_noise(ii);%����1*3000����������Ҳ�������
    end

    w_lms = randn(m, 1);
    w_lms0 = w_lms; 
    w_lms1 = w_lms;
    w_lms2 = w_lms;
    w_lms3 = w_lms;
    w_lms4 = w_lms;
    w_lms5 = w_lms;
    w_lms6 = w_lms;
    w_lms7 = w_lms;
    w_lms8 = w_lms;
    w_lms9 = w_lms;
    w_lms10 = w_lms;
    w_lms11 = w_lms;
    w_GC_LMS = w_lms;
    en = zeros(L,1);%L=10��
    
     %% GMCC_LMS
    alpha_c = 3;
%     lambda = 1/(beta^alpha);
     lambda_c = 0.015;
    
    for ii = 1 : LEN
        Err_GMCC_LMS(mm,ii) = (wo- w_GC_LMS)' *  (wo - w_GC_LMS);
        dn = DD(ii);
        un = UU(:,ii);
        en = dn - w_GC_LMS' * un;
        exp_mid = exp(-1*lambda_c*(abs(en))^alpha_c);
        w_GC_LMS = w_GC_LMS + mu_gmcc*exp_mid*(abs(en))^(alpha_c-1)*sign(en)*un;
%         w_GMCC_LMS = w_GMCC_LMS + mu1 * (en) * un;
        %    w_LMS = w_LMS + mu1 * tanh(en) * un;
    end
    
    %%LMS
    for ii=1:LEN
         err_lms0(mm,ii)=(w_lms0-wo)'*(w_lms0-wo);
         uk=UU(:,ii);
         ek=DD(ii)-w_lms0'*uk;
         w_lms0=w_lms0+mu_lms*ek*uk;
    end
     %%LMM
    for ii=1:LEN
         err_lms11(mm,ii)=(w_lms11-wo)'*(w_lms11-wo);
         uk=UU(:,ii);
         ek=DD(ii)-w_lms11'*uk;
         q = 0;
         a = 3;b = 4;c = 5;
         if ek < a
             q = ek;
         elseif ek >= a && ek < b
             q = a * abs(ek) / ek;
         elseif ek >= b && ek < c
             q = (abs(ek) - c)*a/(b-c)*(abs(ek) / ek);
         elseif ek >= c
             q = 0;
         end
         q = q / ek;
         w_lms11=w_lms11+mu_lmm*ek*uk*q;
    end
    %%MEE
    for ii = L : LEN
        err_lms(mm,ii-L + 1) = (w_lms-wo)'*(w_lms-wo);%mmΪ�����ѭ������1-200�����վ�Ҫ������ͼ��
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_lms' * UU(:,ii - L + jj);
            %en���������ʼ�մ洢10�����1-10��2-11��3-12......
        end 
        for i = 1 : L
            for j = 1 : L
                uij = UU(:,ii - L + i) - UU(:,ii - L + j);
                eij = en(i) - en(j);
                tmp = 1/(L^2*sigma^2) * exp(-(eij^2)/2/sigma^2) * eij * uij;
                w_lms = w_lms + mu_mee * tmp;%mul=0.2;sigma=2;
            end           
        end
    end
    
%% robust5:Hampel ��
    en = zeros(L,1);
    for ii = L : LEN
        err_lms5(mm,ii-L + 1) = (w_lms5-wo)'*(w_lms5-wo);
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_lms5' * UU(:,ii - L + jj);
        end 
        u = 0.6745 * en / median(abs(en - median(en)));
        %median(A)�������Aÿ�е���λ����������Ϊż������Ϊ���м�����ֵ��ƽ��ֵ��
        %en - median(en)��enΪ10*1����������������ÿ��Ԫ�ض���ȥ����λ��,�Ӿ���ֵ��ʹÿ��Ԫ�ض���Ϊ�Ǹ�����
        %u�ķ�ĸ������һ������uΪen�ı�����10*1��������
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
                %���Ȩֵ�ǳ�����en�ϣ�����Ȩ��eij=en(i) - en(j);
                eij = jw5(:,i) - jw5(:,j);
                jw = UU(:,ii - L + i) - UU(:,ii - L + j);%jw����uij��
                tmp =1/(L^2*sigma^2) * exp(-(eij^2)/2/sigma^2) * eij * jw;
                w_lms5 = w_lms5 + mu_hampel * tmp;
            end
        end
    end
 %% robust1:Andrew's��
   en = zeros(L,1);
    for ii = L : LEN
        err_lms1(mm,ii-L + 1) = (w_lms1-wo)'*(w_lms1-wo);
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_lms1' * UU(:,ii - L + jj);
        end 
        u = 0.6745 * en / median(abs(en - median(en)));
        c1= 1.5;
        w1 = zeros(L);
        jw1 = zeros(1,L);
        for jj = 1 : L
            if (u(jj) / c1 == 0)
                w1(jj) = 1;
                jw1(:,jj) = en(jj);
            elseif (pi< abs(u(jj)/c1))
                w1(jj) = 0;
                jw1(:,jj) = 0;
            elseif (abs(u(jj)/c1) <=pi)
                w1(jj) = sin(abs(u(jj)))/abs(u(jj));
                jw1(:,jj) = en(jj) * w1(jj);
            end
        end
        for i = 1 : L
            for j = 1 : L   
                %���Ȩֵ�ǳ�����en�ϣ�����Ȩ��eij=en(i) - en(j);
                eij = jw1(:,i) - jw1(:,j);
                jw = UU(:,ii - L + i) - UU(:,ii - L + j);%jw����uij��
                tmp = 1/(L^2*sigma^2) * exp(-(eij^2)/2/sigma^2) * eij * jw;
                w_lms1 = w_lms1 + mu_andrews * tmp;
            end
        end
    end
    %% robust2:Biweight ��
%    en = zeros(L,1);
%     for ii = L : LEN
%         err_lms2(mm,ii-L + 1) = (w_lms2-wo)'*(w_lms2-wo);
%         for jj = 1 : L
%             en(jj) = DD(ii - L + jj) - w_lms2' * UU(:,ii - L + jj);
%         end 
%         u = 0.6745 * en / median(abs(en - median(en)));
%         c2= 4.685;
%         w2 = zeros(L);
%         jw2 = zeros(1,L);
%         for jj = 1 : L
%             if (abs(u(jj)/c2)>1)
%                 w2(jj) = 0;
%                 jw2(:,jj) = 0;
%             elseif (abs(u(jj)/c2) <=1)
%                 w2(jj) =1-(abs(u(jj)/c2)^4);
%                 jw2(:,jj) = en(jj) * w2(jj);
%             end
%         end
%         for i = 1 : L
%             for j = 1 : L   
%                 %���Ȩֵ�ǳ�����en�ϣ�����Ȩ��eij=en(i) - en(j);
%                 eij = jw2(:,i) - jw2(:,j);
%                 jw = UU(:,ii - L + i) - UU(:,ii - L + j);%jw����uij��
%                 tmp = 1/L^2 * exp(-(eij^2)/2/sigma^2) * eij * jw;
%                 w_lms2 = w_lms2 + mu1 * tmp;
%             end
%         end
%     end
%      %% robust3:Cauchy ��
%    en = zeros(L,1);
%     for ii = L : LEN
%         err_lms3(mm,ii-L + 1) = (w_lms3-wo)'*(w_lms3-wo);
%         for jj = 1 : L
%             en(jj) = DD(ii - L + jj) - w_lms3' * UU(:,ii - L + jj);
%         end 
%         u = 0.6745 * en / median(abs(en - median(en)));
%         c3= 2.385;
%         w3 = zeros(L);
%         jw3 = zeros(1,L);
%         for jj = 1 : L
%                 w3(jj) =1/(1+u(jj)/c3^2);
%                 jw3(:,jj) = en(jj) * w3(jj);
%         end
%         for i = 1 : L
%             for j = 1 : L   
%                 %���Ȩֵ�ǳ�����en�ϣ�����Ȩ��eij=en(i) - en(j);
%                 eij = jw3(:,i) - jw3(:,j);
%                 jw = UU(:,ii - L + i) - UU(:,ii - L + j);%jw����uij��
%                 tmp = 1/L^2 * exp(-(eij^2)/2/sigma^2) * eij * jw;
%                 w_lms3 = w_lms3 + mu1 * tmp;
%             end
%         end
%     end
     %% robust6:Huber��
   en = zeros(L,1);
    for ii = L : LEN
        err_lms6(mm,ii-L + 1) = (w_lms6-wo)'*(w_lms6-wo);
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_lms6' * UU(:,ii - L + jj);
        end 
        u = 0.6745 * en / median(abs(en - median(en)));
        c6= 1.5;
        w6 = zeros(L);
        jw6 = zeros(1,L);
        for jj = 1 : L
            if (abs(u(jj))>c6)
                w6(jj) = c6/abs(u(jj));
                jw6(:,jj) = en(jj) * w6(jj);
            elseif (abs(u(jj))<=c6)
                w6(jj) =1;
                jw6(:,jj) = en(jj) * w6(jj);
            end
        end
        for i = 1 : L
            for j = 1 : L   
                %���Ȩֵ�ǳ�����en�ϣ�����Ȩ��eij=en(i) - en(j);
                eij = jw6(:,i) - jw6(:,j);
                jw = UU(:,ii - L + i) - UU(:,ii - L + j);%jw����uij��
                tmp = 1/(L^2*sigma^2) * exp(-(eij^2)/2/sigma^2) * eij * jw;
                w_lms6 = w_lms6 + mu_huber * tmp;
            end
        end
    end
     %% robust7:Logistic��
%    en = zeros(L,1);
%     for ii = L : LEN
%         err_lms7(mm,ii-L + 1) = (w_lms7-wo)'*(w_lms7-wo);
%         for jj = 1 : L
%             en(jj) = DD(ii - L + jj) - w_lms7' * UU(:,ii - L + jj);
%         end 
%         u = 0.6745 * en / median(abs(en - median(en)));
%         c7= 1.205;
%         w7 = zeros(L);
%         jw7 = zeros(1,L);
%         for jj = 1 : L
%             if (abs(u(jj))==0)
%                 w7(jj) = 1;
%                 jw7(:,jj) = en(jj) * w7(jj);
%             elseif (abs(u(jj))~=0)
%                 w7(jj) =tan(abs(u(jj)/c7))/abs(u(jj)/c7);
%                 jw7(:,jj) = en(jj) * w7(jj);
%             end
%         end
%         for i = 1 : L
%             for j = 1 : L   
%                 %���Ȩֵ�ǳ�����en�ϣ�����Ȩ��eij=en(i) - en(j);
%                 eij = jw7(:,i) - jw7(:,j);
%                 jw = UU(:,ii - L + i) - UU(:,ii - L + j);%jw����uij��
%                 tmp = 1/L^2 * exp(-(eij^2)/2/sigma^2) * eij * jw;
%                 w_lms7 = w_lms7 + mu1 * tmp;
%             end
%         end
%     end
     %% robust8:Median��
%    en = zeros(L,1);
%     for ii = L : LEN
%         err_lms8(mm,ii-L + 1) = (w_lms8-wo)'*(w_lms8-wo);
%         for jj = 1 : L
%             en(jj) = DD(ii - L + jj) - w_lms8' * UU(:,ii - L + jj);
%         end 
%         u = 0.6745 * en / median(abs(en - median(en)));
%         w8 = zeros(L);
%         jw8 = zeros(1,L);
%         for jj = 1 : L
%             if (abs(u(jj))~=0)
%                 w8(jj) = 1/abs(u(jj));
%                 jw8(:,jj) = en(jj) * w8(jj);
%             elseif (abs(u(jj))==0)
%                 w8(jj) =1/0.00001;
%                 jw8(:,jj) = en(jj) * w8(jj);
%             end
%         end
%         for i = 1 : L
%             for j = 1 : L   
%                 %���Ȩֵ�ǳ�����en�ϣ�����Ȩ��eij=en(i) - en(j);
%                 eij = jw8(:,i) - jw8(:,j);
%                 jw = UU(:,ii - L + i) - UU(:,ii - L + j);%jw����uij��
%                 tmp = 1/L^2 * exp(-(eij^2)/2/sigma^2) * eij * jw;
%                 w_lms8 = w_lms8 + mu1 * tmp;
%             end
%         end
%     end

     %% robust9:Talworth��
%    en = zeros(L,1);
%     for ii = L : LEN
%         err_lms9(mm,ii-L + 1) = (w_lms9-wo)'*(w_lms9-wo);
%         for jj = 1 : L
%             en(jj) = DD(ii - L + jj) - w_lms9' * UU(:,ii - L + jj);
%         end 
%         u = 0.6745 * en / median(abs(en - median(en)));
%         c9= 2.795;
%         w9 = zeros(L);
%         jw9 = zeros(1,L);
%         for jj = 1 : L
%             if (abs(u(jj))>c9)
%                 w9(jj) =0;
%                 jw9(:,jj) = en(jj) * w9(jj);
%             elseif (abs(u(jj))<=c9)
%                 w9(jj) =1;
%                 jw9(:,jj) = en(jj) * w9(jj);
%             end
%         end
%         for i = 1 : L
%             for j = 1 : L   
%                 %���Ȩֵ�ǳ�����en�ϣ�����Ȩ��eij=en(i) - en(j);
%                 eij = jw9(:,i) - jw9(:,j);
%                 jw = UU(:,ii - L + i) - UU(:,ii - L + j);%jw����uij��
%                 tmp = 1/L^2 * exp(-(eij^2)/2/sigma^2) * eij * jw;
%                 w_lms9 = w_lms9 + mu1 * tmp;
%             end
%         end
%     end

     %% robust10:Welsch��
%    en = zeros(L,1);
%     for ii = L : LEN
%         err_lms10(mm,ii-L + 1) = (w_lms10-wo)'*(w_lms10-wo);
%         for jj = 1 : L
%             en(jj) = DD(ii - L + jj) - w_lms10' * UU(:,ii - L + jj);
%         end 
%         u = 0.6745 * en / median(abs(en - median(en)));
%         c10= 2.985;
%         w10 = zeros(L);
%         jw10 = zeros(1,L);
%         for jj = 1 : L
%                 w10(jj) =exp(-0.5*u(jj)/c10^2);
%                 jw10(:,jj) = en(jj) * w10(jj);
%         end
%         for i = 1 : L
%             for j = 1 : L   
%                 %���Ȩֵ�ǳ�����en�ϣ�����Ȩ��eij=en(i) - en(j);
%                 eij = jw10(:,i) - jw10(:,j);
%                 jw = UU(:,ii - L + i) - UU(:,ii - L + j);%jw����uij��
%                 tmp = 1/L^2 * exp(-(eij^2)/2/sigma^2) * eij * jw;
%                 w_lms10 = w_lms10 + mu1 * tmp;
%             end
%         end
%     end
     disp(mm)
end
toc


figure(1),hold on;
box on;
wid = 2;
MarkerSize = 8;
Indices = 50;

% ���޵���ɫ����
color_lms0 = [0.85 0.33 0.10]; % ��ɫ
color_lms11 = [0.93 0.69 0.13]; % ��ɫ
color_gmcc = [0.49 0.18 0.56]; % ��ɫ
color_lms = [0.6350 0.0780 0.1840]; % ��ɫ

plot(10*log10(mean(err_lms0)), '-', 'color', color_lms0, 'LineWidth', wid, 'Marker', 'o', 'MarkerSize', MarkerSize, 'MarkerFaceColor', color_lms0, 'MarkerIndices', 1:Indices:LEN)
plot(10*log10(mean(err_lms11)), '-', 'color', color_lms11, 'LineWidth', wid, 'Marker', 's', 'MarkerSize', MarkerSize, 'MarkerFaceColor', color_lms11, 'MarkerIndices', 1:Indices:LEN)
plot(10*log10(mean(Err_GMCC_LMS)), '-', 'color', color_gmcc, 'LineWidth', wid, 'Marker', '^', 'MarkerSize', MarkerSize, 'MarkerFaceColor', color_gmcc, 'MarkerIndices', 1:Indices:LEN)
plot(10*log10(mean(err_lms)), '-c', 'LineWidth', wid, 'Marker', 'd', 'MarkerSize',MarkerSize , 'MarkerFaceColor', 'c', 'MarkerIndices', 1:Indices:LEN)
plot(10*log10(mean(err_lms1)), '-r', 'LineWidth', wid, 'Marker', 'p', 'MarkerSize',MarkerSize , 'MarkerFaceColor', 'r', 'MarkerIndices', 1:Indices:LEN)
plot(10*log10(mean(err_lms6)), '-g', 'LineWidth', wid, 'Marker', 'h', 'MarkerSize', MarkerSize, 'MarkerFaceColor', 'g', 'MarkerIndices', 1:Indices:LEN)
plot(10*log10(mean(err_lms5)), '-b', 'LineWidth', wid, 'Marker', '*', 'MarkerSize', MarkerSize, 'MarkerFaceColor', 'b', 'MarkerIndices', 1:Indices:LEN)


h=legend(['LMS ({\eta\rm}=' num2str(mu_lms) ') '],['LMM (\eta=' num2str(mu_lmm) ', \Delta_1=3, \Delta_2=4, \Delta_3=5)']...
,['GMCC ({\eta=' num2str(mu_gmcc) '}, \alpha=3, \lambda=0.015)'],['MEE ({\eta=' num2str(mu_mee) '}, L=15, \sigma=' num2str(sigma) ')']...
,['MMEE-Andrews ({\eta=' num2str(mu_andrews) '}, L=15, \sigma=' num2str(sigma) ', \Delta_1=1.5)'],['MMEE-Huber ({\eta=' num2str(mu_huber) '}, L=15, \sigma=' num2str(sigma) ', \Delta_1=1.5)']...
,['MMEE-Hample ({\eta=' num2str(mu_hampel) '}, L=15, \sigma=' num2str(sigma) ', \Delta_1=0.5, \Delta_2=2, \Delta_3=4)']);
set(h,'FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',24)
xlabel('Iteration');
ylabel('MSD');


