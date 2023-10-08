function [err_lms1] = F_MEE_andrew(LEN,mu_andrews,wo,w_lms1,DD,UU,sigma,L)
%% robust1:Andrew's
   en = zeros(L,1);
    for ii = L : LEN
        err_lms1(ii-L + 1) = (w_lms1-wo)'*(w_lms1-wo);
        for jj = 1 : L
            en(jj) = DD(ii - L + jj) - w_lms1' * UU(:,ii - L + jj);
        end 
        u = 0.6745 * en / median(abs(en - median(en)));
        c1= 1;
        w1 = zeros(L);
        jw1 = zeros(1,L);
        for jj = 1 : L
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
        for i = 1 : L
            for j = 1 : L   
               
                eij = jw1(:,i) - jw1(:,j);
                jw = UU(:,ii - L + i) - UU(:,ii - L + j);%jw¾ÍÊÇuij£»
%                 tmp = 1/(L^2*sigma^2)  * (1/sqrt(2*pi)/sigma)* exp(-(eij^2)/2/sigma^2) * eij * jw;
                tmp =  exp(-(eij^2)/2/sigma^2) * eij * jw;
                w_lms1 = w_lms1 + mu_andrews * tmp;
            end
        end
    end