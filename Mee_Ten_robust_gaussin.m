%mee -lms
clear all;
close all;
LEN = 3000;
m = 5;
L=20;
sigma1 = 2;
sigma2 = 10;
sigma3 = 20;
sigma4 = 100;
mu0 = 0.045;
mu1 = 0.7;
mu2 = 3;
mu3 = 11;
mu4 = 400;
mu5 = 2;
mu6 = 9;
mu7 = 35;
mu8 = 400;
tic
for mm = 1 : 10
    e_greedy = rand(1, LEN);%���ȷֲ�
    %vv = randn(1,LEN) * 0.1;
    vv=randn(1,LEN)*0.1;
    v1 = randn(1,LEN) * 0.1;
    v2 = randn(1,LEN) * 10;
    VV = (e_greedy <= 0.9).*v1 + (e_greedy >0.9).*v2; %non-gaussian-noise
    VV(1:1000) = v1(1:1000);
    VV(2000:3000) = v1(2000:3000);
    wo1 = randn(m, 1);%��ʼw0Ϊ30*1
    UU = randn(m, LEN);%uuΪ���룬������Ϊ30*1����ȡ3000�����ݣ�
    alpha_noise = alpha_stable_noise(1.1,0.1,0,0,LEN);%alpha=1.1 gamma=0.1
    varNoise = 1;
    alpha_noise = sqrt(varNoise) * alpha_noise;
    for ii = 1 : 1000
        DD(ii) = wo1' * UU(:,ii) + VV(ii);%����1*3000����������Ҳ�������
    end
    wo2 = randn(m, 1);%��ʼw0Ϊ30*1
    for ii = 1001 : 2000
        DD(ii) = wo2' * UU(:,ii) + VV(ii);%����1*3000����������Ҳ�������    
    end
    wo3 = randn(m, 1);%��ʼw0Ϊ30*1
    for ii = 2001 : 3000
        DD(ii) = wo3' * UU(:,ii) + VV(ii);%����1*3000����������Ҳ�������    
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
    % lms
    for ii=L:LEN
         if ii <= 1000
            wo = wo1;
         elseif ii > 1000 && ii <= 2000
            wo = wo2;
         else
            wo = wo3;
        end
         err_lms0(mm,ii-L+1)=(w_lms0-wo)'*(w_lms0-wo);
         uk=UU(:,ii);
         ek=DD(ii)-w_lms0'*uk;
         w_lms0=w_lms0+0.05*ek*uk;
    end
    
    %%MEE
    en = zeros(L,1);%L=10��   
    for ii = L : LEN
        if ii <= 1000
            wo = wo1;
         elseif ii > 1000 && ii <= 2000
            wo = wo2;
         else
            wo = wo3;
        end
        err_mee1(mm,ii-L + 1) = (w_mee1-wo)'*(w_mee1-wo);%mmΪ�����ѭ������1-200�����վ�Ҫ������ͼ��
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_mee1' * UU(:,ii - L + jj);
            %en���������ʼ�մ洢10�����1-10��2-11��3-12......
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
    en = zeros(L,1);%L=10��   
    for ii = L : LEN
        if ii <= 1000
            wo = wo1;
         elseif ii > 1000 && ii <= 2000
            wo = wo2;
         else
            wo = wo3;
        end
        err_mee4(mm,ii-L + 1) = (w_mee4-wo)'*(w_mee4-wo);%mmΪ�����ѭ������1-200�����վ�Ҫ������ͼ��
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_mee4' * UU(:,ii - L + jj);
            %en���������ʼ�մ洢10�����1-10��2-11��3-12......
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
    
    %% robust5:Hampel ��

    en = zeros(L,1);
    for ii = L : LEN
        if ii <= 1000
            wo = wo1;
         elseif ii > 1000 && ii <= 2000
            wo = wo2;
         else
            wo = wo3;
        end
        err_hampel4(mm,ii-L + 1) = (w_hampel4-wo)'*(w_hampel4-wo);
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_hampel4' * UU(:,ii - L + jj);
        end 
        u = 0.6745 * en / median(abs(en - median(en)));
        %median(A)�������Aÿ�е���λ����������Ϊż������Ϊ���м�����ֵ��ƽ��ֵ��
        %en - median(en)��enΪ10*1����������������ÿ��Ԫ�ض���ȥ����λ��,�Ӿ���ֵ��ʹÿ��Ԫ�ض���Ϊ�Ǹ�����
        %u�ķ�ĸ������һ������uΪen�ı�����10*1��������
        aa = 2;
        bb = 4;
        cc = 8;
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

plot(10*log10(mean(err_lms0)),'-b','LineWidth',wid)
% plot(10*log10(mean(err_mee1)),'-.','color',[0.7451 0.7451 0.7451],'LineWidth',wid)
plot(10*log10(mean(err_mee1)),'-g','LineWidth',wid)
% plot(10*log10(mean(err_mee3)),'-.b','LineWidth',wid)
plot(10*log10(mean(err_mee4)),'-c','LineWidth',wid)
% plot(10*log10(mean(err_hampel1)),'-.r','LineWidth',wid)
% plot(10*log10(mean(err_hampel2)),'-.k','LineWidth',wid)
% plot(10*log10(mean(err_hampel3)),'-.m','LineWidth',wid)
plot(10*log10(mean(err_hampel4)),'-r','LineWidth',wid)
%h=legend('LMS','MEE1 (\eta=1.1, L=15, {\bf\sigma=5})','MEE2 ({\eta=3}, L=15, {\bf\sigma=10})','MEE3 ({\eta=11}, L=15, {\bf\sigma=20})','MEE4 ({\eta=600}, L=15, {\bf\sigma=100})','MMEE-Hample1 ({\eta=2},  L=15, {\bf\sigma=5}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4)','MMEE-Hample2 ({\eta=9}, L=15, {\bf\sigma=10}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4)','MMEE-Hample3 ({\eta=35}, L=15, {\bf\sigma=20}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4)','MMEE-Hample4 ({\eta=1000}, L=15, {\bf\sigma=100}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4)');
h=legend('LMS (\eta=0.045)','MEE (\eta=0.7, L=15, {\sigma=2})','MEE ({\eta=400}, L=15, {\sigma=100})','MMEE-Hample ({\eta=400}, L=15, {\sigma=100}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4)');
set(h,'FontName','Times New Roman','FontSize',15,'FontWeight','normal');
xlabel('Iteration');
ylabel('MSD');

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
% legend('MEE-normal','1Andrew��s ','2Biweight','3Cauchy','5Hampel');
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

