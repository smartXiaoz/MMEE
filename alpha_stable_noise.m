function SaS = alpha_stable_noise(alpha,gam,beta,delta,N_train)
    V = -pi/2 + pi.*rand(N_train,1);
    W = exprnd(1,N_train,1);
    if alpha ~= 1
        S_ab = (1+beta^2*tan(alpha*pi/2)^2)^(1/(2*alpha));
        B_ab = atan(beta*tan(alpha*pi/2))/alpha;
        X = S_ab.*sin(alpha.*(V + B_ab))./(cos(V).^(1/alpha)).*( cos(V-alpha.*(V + B_ab))./W  ).^((1-alpha)/alpha);
        SaS = gam.*X + delta;
    else
        X = 2/pi*( (pi/2 + beta.*V).*tan(V) - beta.*log(pi/2.*W.*cos(V)./(pi/2 + beta*V))  );
        SaS = gam.*X + 2/pi*beta*gam + delta;
    end
end