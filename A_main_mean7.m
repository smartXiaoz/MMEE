%mee -lms
clear all;
close all;
LEN = 5000;

L=10;
sigma1 = 2.5;
sigma2 = 3; 
% sigma3 = 3.5;
% sigma4 = 4;
mu5 = 0.08;
mu6 = 0.02;
mu7 = 0.05;
mu8 = 0.3;

[err_hampel1,stdeee1] = F_TH_MMEE(LEN,sigma1,mu5,10,100);
[err_hampel2,stdeee2] = F_TH_MMEE(LEN,sigma2,mu6,10,100);

[err_hampel3,stdeee3] = F_TH_MMEE(LEN,sigma2,mu7,15,100);
[err_hampel4,stdeee4] = F_TH_MMEE(LEN,3.5,mu8,5,100);

figure,hold on;
box on;
wid = 2;
plot(10*log10(mean(err_hampel1)),'-r','LineWidth',wid)
plot(10*log10(mean(err_hampel2)),'-k','LineWidth',wid)
plot(10*log10(mean(err_hampel3)),'-g','LineWidth',wid)
plot(10*log10(mean(err_hampel4)),'-b','LineWidth',wid)
wid = 2.5;
plot(10*log10(ones(1,LEN)*mean(stdeee1)),'-.r','LineWidth',wid)
plot(10*log10(ones(1,LEN)*mean(stdeee2)),'-.k','LineWidth',wid)
plot(10*log10(ones(1,LEN)*mean(stdeee3)),'-.g','LineWidth',wid)
plot(10*log10(ones(1,LEN)*mean(stdeee4)),'-.b','LineWidth',wid)

h1 = 'MMEE-Hample1 ({\eta=0.08},  L=10, {\sigma=2.5})';
h2 = 'MMEE-Hample2 ({\eta=0.02},  L=10, {\sigma=3})';
h3 = 'MMEE-Hample3 ({\eta=0.05},  L=15, {\sigma=3})';
h4 = 'MMEE-Hample4 ({\eta=0.3},  L=5, {\sigma=3.5})';
h=legend(h1,h2,h3,h4,'TH-MMEE-Hample1','TH-MMEE-Hample2','TH-MMEE-Hample3','TH-MMEE-Hample4');
set(h,'FontName','Times New Roman');
xlabel('Iteration');
ylabel('MSD');

TH_1 = 10*log10(mean(stdeee1))
TH_2 = 10*log10(mean(stdeee2))
TH_3 = 10*log10(mean(stdeee3))
TH_4 = 10*log10(mean(stdeee4))
E_Mm1 = mean(10* log10(mean(err_hampel1(end-1000:end))))
E_Mm2 = mean(10* log10(mean(err_hampel2(end-1000:end))))
E_Mm3 = mean(10* log10(mean(err_hampel3(end-1000:end))))
E_Mm4 = mean(10* log10(mean(err_hampel4(end-1000:end))))
