function [err_lms5] = F_MEE_hampel(LEN,mu_hampel,wo,w_lms5,DD,UU,sigma,L)
%% robust5:Hampel 
    en = zeros(L,1);
    for ii = L : LEN
        err_lms5(ii-L + 1) = (w_lms5-wo)'*(w_lms5-wo);
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_lms5' * UU(:,ii - L + jj);
        end 
        u = 0.6745 * en / median(abs(en - median(en)));
        %median(A)�������Aÿ�е���λ����������Ϊż������Ϊ���м�����ֵ��ƽ��ֵ��
        %en - median(en)��enΪ10*1����������������ÿ��Ԫ�ض���ȥ����λ��,�Ӿ���ֵ��ʹÿ��Ԫ�ض���Ϊ�Ǹ�����
        %u�ķ�ĸ������һ������uΪen�ı�����10*1��������
        aa = 0.5;
        bb = 2;
        cc = 4;
        w5 = zeros(L);
        jw5 = zeros(1,L);
        for jj = 1 : L
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
        for i = 1 : L
            for j = 1 : L                 
                eij = jw5(:,i) - jw5(:,j);
                jw = UU(:,ii - L + i) - UU(:,ii - L + j);%jw����uij��
%                 tmp =1/(L^2*sigma^2)  * (1/sqrt(2*pi)/sigma)* exp(-(eij^2)/2/sigma^2) * eij * jw;
                tmp =  exp(-(eij^2)/2/sigma^2) * eij * jw;
                w_lms5 = w_lms5 + mu_hampel * tmp;
            end
        end
    end