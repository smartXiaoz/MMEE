clear all;
close all;

clc

LEN = 6000;
mm = 1;
noise = 0.1;
for pr = [0.75 0.85 0.999]
vv = randn(1,LEN) * noise;
v1=randn(1,LEN)*noise; v2=randn(1,LEN)*10;
rp=rand(1,LEN);   
vv = (rp<=pr).*v1 + (rp>pr).*v2;

MEE = zeros(1, LEN);
M_QIP = zeros(1, LEN);
M_QIP2 = zeros(1, LEN);
L = 20;
sigma = 5;
for ii = L : LEN
    for k = 1 : L
        en(k) = vv(ii - k + 1);
    end
    for i = 1 : L
        for j = 1 : L
            x = en(i) - en(j);
            MEE(ii) = MEE(ii) + 1/(sqrt(2*pi)*sigma)*exp(-x^2/2/sigma^2);
        end
    end
    EQIP(mm,ii-L+1) = MEE(ii) / L^2;
%     E1(mm,ii-L+1) = MEE(ii) / L^2;
end


aa = 2;
bb = 4;
cc = 8;
for ii = L : LEN
    for k = 1 : L
        en(k) = vv(ii - k + 1);
    end
    
    u = 0.6745 * en / median(abs(en - median(en))); 
    for jj = 1 : L
        if abs(u(jj)) <= aa
            w5(jj) = 1; 
        elseif (aa < abs(u(jj))) && (abs(u(jj)) <=bb)
            w5(jj) = aa / abs(u(jj));   
            en(jj) = en(jj) * w5(jj);            
        elseif (bb < abs(u(jj))) && (abs(u(jj)) <=cc)
            w5(jj) = (aa / abs(u(jj))) * ((cc-abs(u(jj)))/(cc-bb));
            en(jj) = en(jj) * w5(jj);            
        elseif abs(u(jj)) > cc
            w5(jj) = 0.000;
            en(jj) = en(jj) * w5(jj);            
        end
    end    
    for i = 1 : L
        for j = 1 : L
            x = en(i) - en(j);
            M_QIP(ii) = M_QIP(ii) + 1/(sqrt(2*pi)*sigma)*exp(-x^2/2/sigma^2);
        end
    end   
    EM_QIP(mm,ii-L+1) = M_QIP(ii) / L^2;
end
mm = mm + 1;
end


% figure(1)
% surf(E1)
% figure(2)
% surf(E2)

figure
hold on
box on
wid = 1;
plot(EM_QIP(1,:),'r','LineWidth',4)
plot(EQIP(1,:),'b','LineWidth',wid)
h=legend('MMEE(\itPr=0.75)','Conventional MEE(\itPr=0.75)');
set(h,'FontName','Times New Roman','FontSize',24,'FontWeight','normal');
% set(gca,'fontsize',24);
ylim([0.04,0.09])
xlabel('Iterations','FontSize',24);
ylabel('Estimate value','FontSize',24);
set(gca,'FontName','Times New Roman');

figure
hold on
box on
wid = 1;
plot(EM_QIP(2,:),'r','LineWidth',4)
plot(EQIP(2,:),'b','LineWidth',wid)
h=legend('MMEE(\itPr=0.85)','Conventional MEE(\itPr=0.85)');
set(h,'FontName','Times New Roman','FontSize',24,'FontWeight','normal');
% set(gca,'fontsize',24);
ylim([0.04,0.09])
xlabel('Iterations','FontSize',24);
ylabel('Estimate value','FontSize',24);
set(gca,'FontName','Times New Roman');
figure
hold on
box on
wid = 1;
plot(EM_QIP(3,:),'r','LineWidth',4)
plot(EQIP(3,:),'b','LineWidth',wid)
h=legend('MMEE(\itPr=0.999)','Conventional MEE(\itPr=0.999)');
set(h,'FontName','Times New Roman','FontSize',24,'FontWeight','normal');
% set(gca,'fontsize',24);
ylim([0.04,0.09])
xlabel('Iterations','FontSize',24);
ylabel('Estimate value','FontSize',24);
