% clear all;
% close all;

% LEN = 10000;
% e_greedy = rand(1, LEN);%均匀分布
% v1 = randn(1,LEN) * 0.1;
% v2 = randn(1,LEN) * 10;
% VV = (e_greedy <= 0.9).*v1 + (e_greedy >0.9).*v2; %non-gaussian-noise
% 
% sigma = 2;
% x = [-10:0.01:10];
% y = exp(-x.^2/(2*sigma^2));

% figure 
% hold on
% plot(x,y)



close all;clear all;
% X = normrnd(0,1,1,10000);%从正态分布中产生10000个均值为0，方差为1 的样本
LEN = 10000;
% X = randn(1,10000)*0.1;
e_greedy = rand(1, LEN);%均匀分布  
stdd = 0.1;
v1 = randn(1,LEN) * stdd;
v2 = randn(1,LEN) * 100;
rp = 0.90;
VV = (e_greedy <= rp).*v1 + (e_greedy > rp).*v2; %non-gaussian-noise
X = VV;
f = -1:0.001:1;%确定横坐标范围



en = X;
u = 0.6745 * en / median(abs(en - median(en)));
aa = 2;
bb = 4;
cc = 8;
w5 = zeros(LEN);
jw5 = zeros(1,LEN);
for jj = 1 : LEN
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
%% Andrew's
c1= 1.5;
w1 = zeros(LEN);
jw1 = zeros(1,LEN);
for jj = 1 : LEN
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

%% huber
c6= 2;
w6 = zeros(LEN);
jw6 = zeros(1,LEN);
for jj = 1 : LEN
    if (abs(u(jj))>c6)
        w6(jj) = c6/abs(u(jj));
        jw6(:,jj) = en(jj) * w6(jj);
    elseif (abs(u(jj))<=c6)
        w6(jj) =1;
        jw6(:,jj) = en(jj) * w6(jj);
    end
end
% index = 1;
% for N = [10 100 1000]
%     for h_1 = [0.05 1 5]
%         p = Parzen(X,h_1,N,f);
%         subplot(3,3,index);
%         plot(f,p,'b'); hold on;
%         legend(['h = ', num2str(h_1), ', N = ',num2str(N)]);
%         plot(f,normpdf(f,0,stdd),'r-');
%         index = index + 1;
%     end
% end 

index = 1;
for N = [100 1000]
    for h_1 = [0.04 0.06]
        p0 = Parzen(X,  h_1,N,f);
        p1 = Parzen(jw5,h_1,N,f);% HAM
        p2 = Parzen(jw6,h_1,N,f);%HUBER
        p3 = Parzen(jw1,h_1,N,f);%AN
        
        subplot(2,2,index);
        hold on;
        wid = 2;
        plot(f,normpdf(f,0,stdd),'-.r','LineWidth',wid);
        plot(f,p0,'LineWidth',wid); 
        plot(f,p1,'LineWidth',wid); 
        plot(f,p2,'LineWidth',wid); 
        plot(f,p3,'LineWidth',wid); 
        h = legend('Normal','MEE','Hampel','Huber','Andrews');
        set(h,'FontName','Times New Roman','FontSize',25,'FontWeight','normal');
        xlabel('e','FontSize',30);
        ylabel('pdf','FontSize',30);
        
%         legend(['h = ', num2str(h_1), ', N = ',num2str(N)]);
%         plot(f,normpdf(f,0,stdd),'r-');
        index = index + 1;
    end
end 
p0 = Parzen(X,0.05,50,f);
p1 = Parzen(jw5,0.05,50,f);% HAM
p2 = Parzen(jw6,0.05,50,f);%HUBER
p3 = Parzen(jw1,0.05,50,f);%AN

figure
hold on
wid = 6;
plot(f,normpdf(f,0,stdd),'-.r','LineWidth',wid);
plot(f,p0,'LineWidth',wid); 
plot(f,p1,'LineWidth',wid); 
plot(f,p2,'LineWidth',wid); 
plot(f,p3,'LineWidth',wid); 

h = legend('Normal','MEE','Hampel','Huber','Andrews');
set(h,'FontName','Times New Roman');
% ylim([-1,1])
xlim([-0.5,0.5])

function p = Parzen(X,h_1,N,x)
% X 是所有的样本
% h 是Parzen窗的窗口大小
% N 是采样的总样本大小
% x 是密度估计的点
% 采用高斯窗口大小

p = zeros(length(x),1); 
h = h_1; 
% N = length(X);
for i = 1:length(x)
    b = 0;
    for j = 1:N
%         ii = randi([1,length(X)]);
        u = (x(i) - X(j))/h;
%         u = (x(i) - X(j));
%         if abs(u) < h
%             b = b + 1;
%         end
        b = b + exp(-u.^2/2)/(sqrt(2*pi));
    end
    p(i) = b/N/h;
end
end