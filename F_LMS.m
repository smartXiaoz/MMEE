function [err_lms0] = F_LMS(LEN,mu_lms,wo,w_lms0,DD,UU)

%%LMS
for ii=1:LEN
     err_lms0(ii)=(w_lms0-wo)'*(w_lms0-wo);
     uk=UU(:,ii);
     ek=DD(ii)-w_lms0'*uk;
     w_lms0=w_lms0+mu_lms*ek*uk;
end