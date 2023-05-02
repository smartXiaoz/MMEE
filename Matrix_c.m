%����MCC������
clear all;
close all;
LEN = 1500;%���г���
m = 5;%ά��
L = 15;%�˲�������������
sigma = 2;%��˹�˺�����bandwidth
mu1 = 0.05;%����
for mm = 1 : 200
    e_greedy = rand(1, LEN);%���ȷֲ�uniform
    %vv = randn(1,LEN) * 0.1;
    v1 = randn(1,LEN) * 0.1;%��˹����������LEN��
    v2 = randn(1,LEN) * 10;%��˹����������LEN��
    VV = (e_greedy <= 0.9).*v1 + (e_greedy >0.9).*v2; %non-gaussian-noise
    %��˹���ģ��
    %�����ˣ���ӦԪ�����
    wo = randn(m, 1);%��ʼ��ʵֵ
    UU = randn(m, LEN);
    w_lms=randn(m,1);%��֤������ֵҪ��ͬ
    w_lms2=w_lms;
    w_lms3=w_lms;
       w_lms4=w_lms;
           w_lms5=w_lms;
    for ii = 1 : LEN
        DD(ii) = wo' * (UU(:,ii) + VV(ii));%DD����1*LEN��Ϊ��֪���������
    end
    
 
    
%number1    
   en = zeros(L,1);%�����
   
    for ii = L : LEN
        err_MEE(mm,ii-L+1) = (w_lms-wo)'*(w_lms-wo);
        for jj = 1:L
            ul(:,jj) = UU(:,ii-jj+1);
            en(jj) = DD(ii-jj+1) - w_lms' * UU(:,ii-jj+1);
        end
    Q = zeros(L, L);
    P = zeros(L, L);
        for i = 1 : L
            for j = 1: L
                x = en(i) - en(j);
                g = exp(-(x^2 )/ 2 / sigma^2);
                Q(i, j) = g;
                P(i, i) = P(i, i) + g;
            end
        end
       
        tmp = (2/ L^2) * ul*(P-Q)*en; 
        w_lms = w_lms + 0.079 * tmp;       
    end
    
    %number 2 a=0.5 b=2 c=2
    en = zeros(L,1);
    for ii = L : LEN
        err_Robust1(mm,ii-L+1) = (w_lms2-wo)'*(w_lms2-wo);
        for jj = 1:L
            ull(:,jj) = UU(:,ii-jj+1);
            en(jj) = DD(ii-jj+1) - w_lms2' * UU(:,ii-jj+1);
        end
        
        
        uu = 0.6745 * en / median(abs(en - median(en)));%10*1 
        aa = 0.5;
        bb = 2;
        cc = 2;
        w7 = zeros(L);
  
        for jj = 1 : L
            if abs(uu(jj)) <= aa
                w7(jj,jj) = 1;        
            elseif (aa < abs(uu(jj))) && (abs(uu(jj)) <=bb)
                w7(jj,jj) = aa / abs(uu(jj));             
            elseif (bb < abs(uu(jj))) && (abs(uu(jj)) <=cc)
                w7(jj,jj) = (aa / abs(uu(jj))) * ((cc-abs(uu(jj)))/(cc-bb));            
            elseif abs(uu(jj)) > cc
                w7(jj,jj) = 0;            
            end
        end
    Q = zeros(L, L);
    P = zeros(L, L);
        for i = 1 : L
            for j = 1: L
                x = w7(i,i)*en(i) - w7(j,j)*en(j);
                g = exp(-(x^2 )/ 2 / sigma^2);
                Q(i, j) = g;
                P(i, i) = P(i, i) + g;
            end
        end
       
        tmp = (2/ L^2)*ull*w7*(P-Q)*w7*en; 
        w_lms2 = w_lms2 + 0.1* tmp;       
    end
    
    %number 3 a=0.5 b=2 c=4
    en = zeros(L,1);
    for ii = L : LEN
        err_Robust2(mm,ii-L+1) = (w_lms3-wo)'*(w_lms3-wo);
        for jj = 1:L
            ull(:,jj) = UU(:,ii-jj+1);
            en(jj) = DD(ii-jj+1) - w_lms3' * UU(:,ii-jj+1);
        end
        
        
        uu = 0.6745 * en / median(abs(en - median(en)));%10*1 
        aa = 0.5;
        bb = 2;
        cc = 4;
        w7 = zeros(L);
  
        for jj = 1 : L
            if abs(uu(jj)) <= aa
                w7(jj,jj) = 1;        
            elseif (aa < abs(uu(jj))) && (abs(uu(jj)) <=bb)
                w7(jj,jj) = aa / abs(uu(jj));             
            elseif (bb < abs(uu(jj))) && (abs(uu(jj)) <=cc)
                w7(jj,jj) = (aa / abs(uu(jj))) * ((cc-abs(uu(jj)))/(cc-bb));            
            elseif abs(uu(jj)) > cc
                w7(jj,jj) = 0;            
            end
        end
    Q = zeros(L, L);
    P = zeros(L, L);
        for i = 1 : L
            for j = 1: L
                x = w7(i,i)*en(i) - w7(j,j)*en(j);
                g = exp(-(x^2 )/ 2 / sigma^2);
                Q(i, j) = g;
                P(i, i) = P(i, i) + g;
            end
        end
       
        tmp = (2/ L^2)*ull*w7*(P-Q)*w7*en; 
        w_lms3 = w_lms3 + 0.09 * tmp;       
    end
    
        %number 4 a=0.5 b=2 c=8
    en = zeros(L,1);
    for ii = L : LEN
        err_Robust3(mm,ii-L+1) = (w_lms4-wo)'*(w_lms4-wo);
        for jj = 1:L
            ull(:,jj) = UU(:,ii-jj+1);
            en(jj) = DD(ii-jj+1) - w_lms4' * UU(:,ii-jj+1);
        end
        
        
        uu = 0.6745 * en / median(abs(en - median(en)));%10*1 
        aa = 0.5;
        bb = 2;
        cc = 8;
        w7 = zeros(L);
  
        for jj = 1 : L
            if abs(uu(jj)) <= aa
                w7(jj,jj) = 1;        
            elseif (aa < abs(uu(jj))) && (abs(uu(jj)) <=bb)
                w7(jj,jj) = aa / abs(uu(jj));             
            elseif (bb < abs(uu(jj))) && (abs(uu(jj)) <=cc)
                w7(jj,jj) = (aa / abs(uu(jj))) * ((cc-abs(uu(jj)))/(cc-bb));            
            elseif abs(uu(jj)) > cc
                w7(jj,jj) = 0;            
            end
        end
    Q = zeros(L, L);
    P = zeros(L, L);
        for i = 1 : L
            for j = 1: L
                x = w7(i,i)*en(i) - w7(j,j)*en(j);
                g = exp(-(x^2 )/ 2 / sigma^2);
                Q(i, j) = g;
                P(i, i) = P(i, i) + g;
            end
        end
       
        tmp = (2/ L^2)*ull*w7*(P-Q)*w7*en; 
        w_lms4 = w_lms4 + 0.098 * tmp;       
    end
            %number 5 a=0.5 b=2 c=10
    en = zeros(L,1);
    for ii = L : LEN
        err_Robust4(mm,ii-L+1) = (w_lms5-wo)'*(w_lms5-wo);
        for jj = 1:L
            ull(:,jj) = UU(:,ii-jj+1);
            en(jj) = DD(ii-jj+1) - w_lms5' * UU(:,ii-jj+1);
        end
        
        
        uu = 0.6745 * en / median(abs(en - median(en)));%10*1 
        aa = 0.5;
        bb = 2;
        cc = 10;
        w7 = zeros(L);
  
        for jj = 1 : L
            if abs(uu(jj)) <= aa
                w7(jj,jj) = 2;        
            elseif (aa < abs(uu(jj))) && (abs(uu(jj)) <=bb)
                w7(jj,jj) = aa / abs(uu(jj));             
            elseif (bb < abs(uu(jj))) && (abs(uu(jj)) <=cc)
                w7(jj,jj) = (aa / abs(uu(jj))) * ((cc-abs(uu(jj)))/(cc-bb));            
            elseif abs(uu(jj)) > cc
                w7(jj,jj) = 0;            
            end
        end
    Q = zeros(L, L);
    P = zeros(L, L);
        for i = 1 : L
            for j = 1: L
                x = w7(i,i)*en(i) - w7(j,j)*en(j);
                g = exp(-(x^2 )/ 2 / sigma^2);
                Q(i, j) = g;
                P(i, i) = P(i, i) + g;
            end
        end
       
        tmp = (2/ L^2)*ull*w7*(P-Q)*w7*en; 
        w_lms5 = w_lms5 + 0.095 * tmp;       
    end
    disp(mm)
end
box on;%box on����ϵ�ұߺ��ϱ��б߿�
figure(1);
% plot(10*log10(mean(err_MCC)));
hold on;
plot(10*log10(mean(err_MEE)));
hold on;
plot(10*log10(mean(err_Robust1)));
hold on;
plot(10*log10(mean(err_Robust2)));
hold on;
plot(10*log10(mean(err_Robust3)));
hold on;
plot(10*log10(mean(err_Robust4)));
legend('MEE','MEE+Hampel(c=2)','MEE+Hampel(c=4)','MEE+Hampel(c=8)','MEE+Hampel(c=10)');
xlabel('iteration');
ylabel('Weight error power');
10*log10(mean(err_MEE(:,1486)))
10*log10(mean(err_Robust1(:,1486)))
10*log10(mean(err_Robust2(:,1486)))
10*log10(mean(err_Robust3(:,1486)))
10*log10(mean(err_Robust4(:,1486)))
% legend('MCC','MEE','MEE+Hampel-Robust');