%mee -lms
clear all;
close all;
LEN = 6000;

L=10;

sigma1 = 3.5; 
% sigma3 = 3.5;
% sigma4 = 4;
emmu = [0.05:0.02:0.5];
mm = 1;
for mu5 = emmu
[err_hampel1,stdeee1] = F_TH_MMEE(LEN,sigma1,mu5,10,20);
f_mean(mm) = mean(mean(err_hampel1(:,end-3000:end)));
t_mean(mm) = mean(stdeee1);
disp(mm)
mm = mm + 1;
end


figure,hold on;
box on;
wid = 6;
plot(emmu,(f_mean),'-r','LineWidth',wid)
plot(emmu,(t_mean),':bO','LineWidth',wid)
h = legend('simulation','theory');
xlabel('\eta');
ylabel('MSD');
