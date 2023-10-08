clear all;
close all;
mu_mee = 0.001;
mu_andrews = 0.001;
mu_huber = 0.001;
mu_hampel = 0.001;
sigma_all = 0.1:0.05:2;
LEN = size(sigma_all, 2);
parfor ii = 1 : LEN
    sigma = sigma_all(ii);
    [Err_mee(ii),Err_mee_hampel(ii),Err_mee_andrew(ii),Err_mee_huber(ii)] = F_sigma(sigma);
    disp(ii)
end

figure 
hold on
box on
wid = 5;
plot(sigma_all,10*log10(Err_mee),'c-s','LineWidth',wid)
plot(sigma_all,10*log10(Err_mee_andrew),'g-o','LineWidth',wid)
plot(sigma_all,10*log10(Err_mee_huber),'r-d','LineWidth',wid)
plot(sigma_all,10*log10(Err_mee_hampel),'b-o','LineWidth',wid)


h=legend(['MEE ({\eta=' num2str(mu_mee) '}, L=15)']...
,['MMEE-Andrews ({\eta=' num2str(mu_andrews) '}, L=15, \Delta_1=1)'],['MMEE-Huber ({\eta=' num2str(mu_huber) '}, L=15, \Delta_1=0.75)']...
,['MMEE-Hample ({\eta=' num2str(mu_hampel) '}, L=15, \Delta_1=0.5, \Delta_2=2, \Delta_3=4)']);
set(h,'FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',24)
xlabel('\sigma','FontSize',24);
ylabel('MSD(dB)','FontSize',24);
