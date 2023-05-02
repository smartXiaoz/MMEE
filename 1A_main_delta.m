%mee -lms
clear all;
close all;
LEN = 1000;
m = 5;
L=15;
sigma1 = 2;
sigma2 = 2;
sigma3 = 2;
sigma4 = 2;

mu1 = 0.45;
mu5 = 0.34;
mu6 = 0.24;
mu7 = 0.29;
mu8 = 0.35;
tic
for mm = 1 : 500
    e_greedy = rand(1, LEN);%���ȷֲ�
    %vv = randn(1,LEN) * 0.1;
    vv=randn(1,LEN)*0.1;
    v1 = randn(1,LEN) * 0.1;
    v2 = randn(1,LEN) * 10;
    VV = (e_greedy <= 0.9).*v1 + (e_greedy >0.9).*v2; %non-gaussian-noise
    
    wo = randn(m, 1);%��ʼw0Ϊ30*1
    UU = randn(m, LEN);%uuΪ���룬������Ϊ30*1����ȡ3000�����ݣ�
    alpha_noise = alpha_stable_noise(1.1,0.1,0,0,LEN);%alpha=1.1 gamma=0.1
    varNoise = 1;
    alpha_noise = sqrt(varNoise) * alpha_noise;
    for ii = 1 : LEN
        DD(ii) = wo' * UU(:,ii) + VV(ii);%����1*3000����������Ҳ�������
%         DD(ii) = wo' * UU(:,ii) +alpha_noise(ii);%����1*3000����������Ҳ�������
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
    
    en = zeros(L,1);%L=10��   
    for ii = L : LEN
        err_mee1(mm,ii-L + 1) = (w_mee1-wo)'*(w_mee1-wo);%mmΪ�����ѭ������1-200�����վ�Ҫ������ͼ��
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_mee1' * UU(:,ii - L + jj);
            ul(:,jj) = UU(:,ii-jj+1);
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
    
    
%% robust5:Hampel ��
    en = zeros(L,1);
    for ii = L : LEN
        err_hampel1(mm,ii-L + 1) = (w_hampel1-wo)'*(w_hampel1-wo);
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_hampel1' * UU(:,ii - L + jj);
        end 
        u = 0.6745 * en / median(abs(en - median(en)));
        %median(A)�������Aÿ�е���λ����������Ϊż������Ϊ���м�����ֵ��ƽ��ֵ��
        %en - median(en)��enΪ10*1����������������ÿ��Ԫ�ض���ȥ����λ��,�Ӿ���ֵ��ʹÿ��Ԫ�ض���Ϊ�Ǹ�����
        %u�ķ�ĸ������һ������uΪen�ı�����10*1��������
        aa = 0.25;
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
                tmp =1/(L^2*sigma1^2) * exp(-(eij^2)/2/sigma1^2) * eij * jw;
                w_hampel1 = w_hampel1 + mu5 * tmp;
            end
        end
    end
    
    
    
    %% robust5:Hampel ��

    en = zeros(L,1);
    for ii = L : LEN
        err_hampel2(mm,ii-L + 1) = (w_hampel2-wo)'*(w_hampel2-wo);
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_hampel2' * UU(:,ii - L + jj);
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
                tmp =1/(L^2*sigma2^2) * exp(-(eij^2)/2/sigma2^2) * eij * jw;
                w_hampel2 = w_hampel2 + mu6 * tmp;
            end
        end
    end
    
    %% robust5:Hampel ��

    en = zeros(L,1);
    for ii = L : LEN
        err_hampel3(mm,ii-L + 1) = (w_hampel3-wo)'*(w_hampel3-wo);
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_hampel3' * UU(:,ii - L + jj);
        end 
        u = 0.6745 * en / median(abs(en - median(en)));
        %median(A)�������Aÿ�е���λ����������Ϊż������Ϊ���м�����ֵ��ƽ��ֵ��
        %en - median(en)��enΪ10*1����������������ÿ��Ԫ�ض���ȥ����λ��,�Ӿ���ֵ��ʹÿ��Ԫ�ض���Ϊ�Ǹ�����
        %u�ķ�ĸ������һ������uΪen�ı�����10*1��������
        aa = 1;
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
                tmp =1/(L^2*sigma3^2) * exp(-(eij^2)/2/sigma3^2) * eij * jw;
                w_hampel3 = w_hampel3 + mu7 * tmp;
            end
        end
    end
    
    %% robust5:Hampel ��

    en = zeros(L,1);
    for ii = L : LEN
        err_hampel4(mm,ii-L + 1) = (w_hampel4-wo)'*(w_hampel4-wo);
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_hampel4' * UU(:,ii - L + jj);
        end 
        u = 0.6745 * en / median(abs(en - median(en)));
        %median(A)�������Aÿ�е���λ����������Ϊż������Ϊ���м�����ֵ��ƽ��ֵ��
        %en - median(en)��enΪ10*1����������������ÿ��Ԫ�ض���ȥ����λ��,�Ӿ���ֵ��ʹÿ��Ԫ�ض���Ϊ�Ǹ�����
        %u�ķ�ĸ������һ������uΪen�ı�����10*1��������
        aa = 1.5;
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
plot(10*log10(mean(err_mee1)),'-b','LineWidth',wid)

plot(10*log10(mean(err_hampel1)),'-r','LineWidth',wid)
plot(10*log10(mean(err_hampel2)),'-k','LineWidth',wid)
plot(10*log10(mean(err_hampel3)),'-m','LineWidth',wid)
plot(10*log10(mean(err_hampel4)),'-','color',[0.7451 0.451 0.451],'LineWidth',wid)
h1 = 'MEE1 (\eta=0.45, L=15, {\sigma=2})';
h2 = 'MMEE-Hampel1 ({\eta=0.34},  L=15, {\sigma=2}, {\bf\Delta_1=0.25}, {\Delta_2=2}, {\Delta_3=4})';
h3 = 'MMEE-Hampel2 ({\eta=0.24},  L=15, {\sigma=2}, {\bf\Delta_1=0.5}, {\Delta_2=2}, {\Delta_3=4})';
h4 = 'MMEE-Hampel3 ({\eta=0.29},  L=15, {\sigma=2}, {\bf\Delta_1=1}, {\Delta_2=2}, {\Delta_3=4})';
h5 = 'MMEE-Hampel4 ({\eta=0.35},  L=15, {\sigma=2}, {\bf\Delta_1=1.5}, {\Delta_2=2}, {\Delta_3=4})';
h=legend(h1,h2,h3,h4,h5);
set(h,'FontName','Times New Roman','FontSize',15,'FontWeight','normal');
xlabel('Iteration');
ylabel('MSD');

E_MEE1 = mean(10* log10(mean(err_mee1(end/2:end))));

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

