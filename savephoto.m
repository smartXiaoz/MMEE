open('../photo/pdf.fig');
labels = get(legend(), 'String'); 
plots = flipud(get(gca, 'children')); 

% Now re-create the legend 
neworder = [5,1,2,3,4]; %ͼ���ı���˳��
box on
legend(plots(neworder), labels(neworder))
