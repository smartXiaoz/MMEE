clear all;
close all;
 

% 生成一些正态分布的随机数
% x = randn(1,10000);
LEN = 10000;
% X = randn(1,10000)*0.1;
e_greedy = rand(1, LEN);%均匀分布   
v1 = randn(1,LEN) * 1;
v2 = randn(1,LEN) * 100;
VV = (e_greedy <= 0.9).*v1 + (e_greedy >0.9).*v2; %non-gaussian-noise
x = VV;
n = LEN;
minx = min(x);
maxx = max(x);
dx = (maxx-minx)/n;
x1 = -3:0.01:3;
n = length(x1);
h=0.5;
f=zeros(1,n);

% N = 10000;
% for j = 1:n
%     for i=1:N
%         f(j)=f(j)+exp(-(x1(j)-x(i))^2/2/h^2)/sqrt(2*pi);
%     end
%     f(j)=f(j)/N/h;
% end
% plot(x1,f);
  
%用系统函数计算比较
[f2,x2] = ksdensity(x);
hold on;
plot(x2,f2,'r'); %红色为参考