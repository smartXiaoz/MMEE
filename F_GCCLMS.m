function [err_GMCC_LMS] = F_GCCLMS(LEN,mu_gmcc,wo,w_GC_LMS,DD,UU)

%% GMCC_LMS
alpha_c = 3;
%     lambda = 1/(beta^alpha);
lambda_c = 0.015;

for ii = 1 : LEN
    err_GMCC_LMS(ii) = (wo- w_GC_LMS)' *  (wo - w_GC_LMS);
    dn = DD(ii);
    un = UU(:,ii);
    en = dn - w_GC_LMS' * un;
    exp_mid = exp(-1*lambda_c*(abs(en))^alpha_c);
    w_GC_LMS = w_GC_LMS + mu_gmcc*exp_mid*(abs(en))^(alpha_c-1)*sign(en)*un;
%         w_GMCC_LMS = w_GMCC_LMS + mu1 * (en) * un;
    %    w_LMS = w_LMS + mu1 * tanh(en) * un;
end