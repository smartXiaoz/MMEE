function [err_lms6] = F_MEE_huber(LEN,mu_huber,wo,w_lms6,DD,UU,sigma,L) 
%% robust6:Huber
   en = zeros(L,1);
    for ii = L : LEN
        err_lms6(ii-L + 1) = (w_lms6-wo)'*(w_lms6-wo);
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_lms6' * UU(:,ii - L + jj);
        end 
        u = 0.6745 * en / median(abs(en - median(en)));
        c6= 0.75;
        w6 = zeros(L);
        jw6 = zeros(1,L);
        for jj = 1 : L
            if (abs(u(jj))>c6)
                w6(jj) = c6/abs(u(jj));
                jw6(:,jj) = en(jj) * w6(jj);
            elseif (abs(u(jj))<=c6)
                w6(jj) =1;
                jw6(:,jj) = en(jj) * w6(jj);
            end
        end
        for i = 1 : L
            for j = 1 : L   
               
                eij = jw6(:,i) - jw6(:,j);
                jw = UU(:,ii - L + i) - UU(:,ii - L + j);%jw¾ÍÊÇuij£»
%                 tmp = 1/(L^2*sigma^2)  * (1/sqrt(2*pi)/sigma)* exp(-(eij^2)/2/sigma^2) * eij * jw;
                tmp =  exp(-(eij^2)/2/sigma^2) * eij * jw;
                w_lms6 = w_lms6 + mu_huber * tmp;
            end
        end
    end