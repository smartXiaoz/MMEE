clear all;
close all;
clc
% 打开 fig 文件
fig = openfig('erle.fig', 'invisible');

% 找到图形对象，根据需求自动调整参数
lines = findobj(fig, 'type', 'line')

% 获取数据
xdata = get(lines, 'XData');
ydata = get(lines, 'YData');
err_mee2 = ydata{6};
err_mee1 = ydata{5};
err_hampel1 = ydata{4};
err_hampel2 = ydata{3};
err_hampel3 = ydata{2};
err_hampel4 = ydata{1};
% % 绘制所有曲线
% figure;
% hold on;
% for i = 1:length(lines)
%     plot(xdata{i}, ydata{i});
% end
close(fig)
LEN = 26800;
figure(1),hold on;
box on;
wid = 2.5;
MarkerSize = 6;
Indices=500;
plot(err_mee2, '-k', 'LineWidth', wid, 'Marker', 'd', 'MarkerSize',MarkerSize , 'MarkerFaceColor', 'k', 'MarkerIndices', 1:Indices:LEN)
plot(err_mee1, '-c', 'LineWidth', wid, 'Marker', 'd', 'MarkerSize',MarkerSize , 'MarkerFaceColor', 'c', 'MarkerIndices', 1:Indices:LEN)
plot(err_hampel1, '-r', 'LineWidth', wid, 'Marker', 'p', 'MarkerSize',MarkerSize , 'MarkerFaceColor', 'r', 'MarkerIndices', 1:Indices:LEN)
plot(err_hampel2, '-g', 'LineWidth', wid, 'Marker', 'h', 'MarkerSize', MarkerSize, 'MarkerFaceColor', 'g', 'MarkerIndices', 1:Indices:LEN)
plot(err_hampel3, '-b', 'LineWidth', wid, 'Marker', '*', 'MarkerSize', MarkerSize, 'MarkerFaceColor', 'b', 'MarkerIndices', 1:Indices:LEN)
plot(err_hampel4, '-m', 'LineWidth', wid, 'Marker', '*', 'MarkerSize', MarkerSize, 'MarkerFaceColor', 'm', 'MarkerIndices', 1:Indices:LEN)

h6 = '(MMEE-Huber)';
h5 = '(MMEE-Andrews)';
h4 = '(MMEE-Hampel)';
h3 = '(MEE)';
h2 = '(GMCC)';
h1 = '(LMS)';
xlim([1,26000])
h=legend(h1,h2,h3,h4,h5,h6);
set(h,'FontName','Times New Roman','FontSize',24,'FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',24,'FontWeight','normal');
xlabel('Iteration');
ylabel('MSD');
