close all;
clear all;
e_greedy = rand(1, 181);%¾ùÔÈ·Ö²¼
v1 = randn(1,181) * 0.55;
    v2 = randn(1,181) * 5;
    VV = (e_greedy <= 0.85).*v1 + (e_greedy >0.85).*v2;
x=1:0.05:10;
y=2*x+1;
y1=2*x+1+1.2*VV;
y2=0.8*x+7;

figure
box on
hold on
axis([0 11.5 0 35])
plot(x,y,'-r','Linewidth',8);
plot(x,y2,'-b','Linewidth',8);
plot(x,y1,'.k','MarkerSize',30);
h = legend('Robust Regression','LS Regression');
set(h,'FontName','Times New Roman','FontSize',25,'FontWeight','normal')
ylim([0,25])
xlabel('x','FontSize',30)
ylabel('y','FontSize',30)
