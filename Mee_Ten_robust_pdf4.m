close all;clear all;
% X = normrnd(0,1,1,10000);%从正态分布中产生10000个均值为0，方差为1 的样本
LEN = 10000;

e_greedy = rand(1, LEN);%均匀分布  
stdd = 0.1;
v1 = randn(2,LEN) * stdd;
v2 = randn(2,LEN) * 100;
VV = (e_greedy <= 0.9).*v1 + (e_greedy >0.9).*v2; %non-gaussian-noise
X = VV;
f = -1:0.01:1;%确定横坐标范围






en = X(1,:);
u = 0.6745 * en / median(abs(en - median(en)));
aa = 1;
bb = 2;
cc = 100;
w5 = zeros(LEN);
jw5 = zeros(2,LEN);
for jj = 1 : LEN
    if abs(u(jj)) <= aa
        w5(jj) = 1;
        jw5(1,jj) = en(jj);
    elseif (aa < abs(u(jj))) && (abs(u(jj)) <=bb)
        w5(jj) = aa / abs(u(jj));
        jw5(1,jj) = en(jj) * w5(jj);
    elseif (bb < abs(u(jj))) && (abs(u(jj)) <=cc)
        w5(jj) = (aa / abs(u(jj))) * ((cc-abs(u(jj)))/(cc-bb));
        jw5(1,jj) = en(jj) * w5(jj);
    elseif abs(u(jj)) > cc
        w5(jj) = 0;
        jw5(1,jj) = 0;
    end
end

en = X(2,:);
u = 0.6745 * en / median(abs(en - median(en)));
w5 = zeros(LEN);
for jj = 1 : LEN
    if abs(u(jj)) <= aa
        w5(jj) = 1;
        jw5(2,jj) = en(jj);
    elseif (aa < abs(u(jj))) && (abs(u(jj)) <=bb)
        w5(jj) = aa / abs(u(jj));
        jw5(2,jj) = en(jj) * w5(jj);
    elseif (bb < abs(u(jj))) && (abs(u(jj)) <=cc)
        w5(jj) = (aa / abs(u(jj))) * ((cc-abs(u(jj)))/(cc-bb));
        jw5(2,jj) = en(jj) * w5(jj);
    elseif abs(u(jj)) > cc
        w5(jj) = 0;
        jw5(2,jj) = 0;
    end
end

h = 0.1;
N = 10;
LEN = length(f);
y = zeros(LEN,LEN);
y1 = zeros(LEN,LEN);
for i = 1 : LEN
    for j = 1 : LEN
        x(1) = f(i);
        x(2) = f(j);      
        for kk = 1 : N
            u = norm((x' - X(:,kk)))/h;
            y(i,j) = y(i,j) + exp(-u.^2/2)/(sqrt(2*pi));
            u1 = norm((x' - jw5(:,kk)))/h;
            y1(i,j) = y1(i,j) + exp(-u1.^2/2)/(sqrt(2*pi));
        end
        y(i,j) = y(i,j)/N/h;
        y1(i,j) = y1(i,j)/N/h;
    end
end




sigma = stdd;
Z = zeros(LEN , LEN);

for row = 1 : 1 : LEN
    for col = 1 : 1 : LEN
        Z( row, col ) = ( f(row)).^2 + ( f(col)).^2; 
    end
end
Z = -Z / ( 2 * sigma^2 ); 
Z = exp(Z) / ( 2 * pi * sigma^2 );


figure; 
hold on 
mesh(f, f, y, 'edgecolor','r','linewidth',0.5)
mesh(f, f, y1, 'edgecolor','b','linewidth',0.5)

figure; 
hold on 
contour(f, f, y, 'edgecolor','r','linewidth',0.5)
contour(f, f, y1, 'edgecolor','b','linewidth',0.5)
figure(2); 
hold on 
mesh(f, f, Z, 'edgecolor','r','linewidth',0.5)




