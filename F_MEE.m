 function [err_lms] = F_MEE(LEN,mu_mee,wo,w_lms,DD,UU,sigma,L)
 
%%MEE
    for ii = L : LEN
        err_lms(ii-L + 1) = (w_lms-wo)'*(w_lms-wo);%mmΪ�����ѭ������1-200�����վ�Ҫ������ͼ��
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_lms' * UU(:,ii - L + jj);
            %en���������ʼ�մ洢10�����1-10��2-11��3-12......
        end 
        for i = 1 : L
            for j = 1 : L
                uij = UU(:,ii - L + i) - UU(:,ii - L + j);
                eij = en(i) - en(j);
%                 tmp = 1/(L^2*sigma^2) * (1/sqrt(2*pi)/sigma)*exp(-(eij^2)/2/sigma^2) * eij * uij;
                tmp = exp(-(eij^2)/2/sigma^2) * eij * uij;
                w_lms = w_lms + mu_mee * tmp;%mul=0.2;sigma=2;
            end           
        end
    end