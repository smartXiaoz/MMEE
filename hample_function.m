clear all;
close all;
a=2;
b=4;
c=8;
u=-10:0.05:10;
A=[];
Ha=[];
Hu=[];
%%1
for uu=u
    if (abs(uu)<=a)
        Ha=[Ha,1];
    elseif (abs(uu)>a&&abs(uu)<=b)
        Ha=[Ha,a/abs(uu)];
    elseif (abs(uu)>b&&abs(uu)<=c)
        Ha=[Ha,(a*((c-abs(uu))/(c-b)))/(abs(uu))];
    elseif (abs(uu)>c)
        Ha=[Ha,0];
    end
    if (uu/1.339==0)
        A=[A,1];
    elseif (abs(uu/1.339)>pi)
        A=[A,0];
    elseif (abs(uu/1.339)<=pi)
        A=[A,(sin(abs(uu/c)))/(abs(uu/c))];
    end
    if (abs(uu)<=1.345)
        Hu=[Hu,1];
    elseif (abs(uu)>1.345)
        Hu=[Hu,1.345/abs(uu)];
    end
end
figure(1);
hold on;
plot(u,A,'r-.','LineWidth',10);
plot(u,Ha,'b-.','LineWidth',10);
plot(u,Hu,'g-.','LineWidth',10);
axis([-10 10 -0 1.3]);
xlabel('x','FontSize',20)
ylabel('y','FontSize',20)
h=legend('Andrews \Delta_1=1.339','Hampel  \Delta_1=2, \Delta_2=4, \Delta_3=8','Huber     \Delta_1=1.345');
set( h, 'FontSize', 22);
%grid on;
box on;