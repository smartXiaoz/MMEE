function [err_lms11] = F_LMM(LEN,mu_lmm,wo,w_lms11,DD,UU)
 %%LMM
for ii=1:LEN
     err_lms11(ii)=(w_lms11-wo)'*(w_lms11-wo);
     uk=UU(:,ii);
     ek=DD(ii)-w_lms11'*uk;
     q = 0;
     a = 3;b = 4;c = 5;
     if ek < a
         q = ek;
     elseif ek >= a && ek < b
         q = a * abs(ek) / ek;
     elseif ek >= b && ek < c
         q = (abs(ek) - c)*a/(b-c)*(abs(ek) / ek);
     elseif ek >= c
         q = 0;
     end
     q = q / ek;
     w_lms11=w_lms11+mu_lmm*ek*uk*q;
end