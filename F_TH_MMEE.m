function [err, th_err] = F_TH_MMEE(LEN,sigma,mu,L,mmm)
for mm = 1 : mmm
m = 5;
e_greedy = rand(1, LEN);
vv=randn(1,LEN)*0.1;
v1 = randn(1,LEN) * 0.1;
v2 = randn(1,LEN) * 10;
VV = (e_greedy <= 0.99).*v1 + (e_greedy >0.99).*v2; %non-gaussian-noise
wo = randn(m, 1);%初始w0为30*1
UU = randn(m, LEN);%uu为输入，理论上为30*1，但取3000组数据；
for ii = 1 : LEN
    DD(ii) = wo' * UU(:,ii) + VV(ii);%产生1*3000的行向量，也即是输出
end
w_hampel1 = randn(m, 1);
%% MVV


%% robust5:Hampel 法
Dm2 = zeros(m,m);
Dm3 = zeros(m,1);
en = zeros(L,1);
DM6 = 0;
DW0 = 0;
DM7 = 0;
DM8 = 0;
DM9 = 0;
DM10 = 0;
for ii = L : LEN
    vvll = zeros(L,1);
    ww = zeros(L,L);
    ul = zeros(m,L);
    err_hampel1(mm,ii-L + 1) = (w_hampel1-wo)'*(w_hampel1-wo);
    for jj = 1 : L
        en(jj) = DD(ii - L + jj) - w_hampel1' * UU(:,ii - L + jj);
        vvll(jj) = VV(:,ii - L + jj);
        ul(:,jj) = UU(:,ii - L + jj);
    end 
    u = 0.6745 * en / median(abs(en - median(en)));
    %median(A)是求矩阵A每列的中位数，若个数为偶数，则为最中间两个值的平均值；
    %en - median(en)：en为10*1的列向量，该向量每个元素都减去其中位数,加绝对值会使每个元素都变为非负数；
    %u的分母最终是一个数，u为en的倍数（10*1列向量）
    aa = 2;
    bb = 4;
    cc = 8;
    w5 = zeros(L,1);
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
    P = zeros(L, L);
    Q = zeros(L, L);
    P1 = zeros(L, L);
    Q1 = zeros(L, L);
    for i = 1 : L
        for j = 1 : L   
            %这个权值是乘在了en上，不加权重eij=en(i) - en(j);
            eij = jw5(:,i) - jw5(:,j);
            eij1 = en(i) - en(j);
            jw = UU(:,ii - L + i) - UU(:,ii - L + j);%jw就是uij；
            g = exp(-(eij^2 )/ 2 / sigma^2);
            g1 = exp(-(eij1^2 )/ 2 / sigma^2);
            Q(i, j) = g;
            P(i, i) = P(i, i) + g;
            Q1(i, j) = g1;
            P1(i, i) = P(i, i) + g1;
            tmp =1/(L^2*sigma^2) * g * eij * jw;
            w_hampel1 = w_hampel1 + mu * tmp;
        end
    end        
    if ii > 3000
        ww = diag(w5);
%         Dm2 = Dm2 + ul*ww*(P-Q)*ww*ul';
%         Dm3 = Dm3 + ul*ww*(P-Q)*ww*vvll;
        DW0 = DW0 + ww;
        DM6 = DM6 + vvll'*ww*(P-Q)*ww*ww*(P-Q)*ww*vvll;
        DM7 = DM7 + ww*(P-Q)*ww*ww*(P-Q)*ww;
        DM8 = DM8 + ww*(P-Q)*ww;
        DM9 = DM9 + vvll'*ww*vvll;
        DM10 = DM10 + (P1-Q1);
    end
end
Mvv = mhampel(VV,LEN);
D6 = 0;
for i = 1 : L-1
    for j = 3000  : LEN 
       vij = Mvv(j) - Mvv(j-i);
       D6 = D6 + exp(-(vij^2 )/ 2 / sigma^2);
    end
end
LEN1 = LEN-3000;
D6 = D6/LEN1;

yita =  1/(L^2*sigma^2);
mu11 = yita*mu;
tmp1 = mean(diag(DW0/LEN1));
eep = DM6/LEN1*m*L*tmp1^2;
% eep = DM6/LEN1*m*L;
% eep = D6^2*m*var(VV)*tmp1^4;
% eep = (mean(diag(DM8/LEN1))^2 + (L-1)*tmp1^2) * var(VV) * L




aaa = mu11^2*eep;       
eeu = D6*tmp1^2*L;    
eeccc = 1-(1-mu11*2*eeu)^2;
stdeeem(mm) = 4*aaa/eeccc;

% eep = ((Dm3'*Dm3)/LEN1);
% yita =  1/(L^2*sigma^2);
% mu11 = yita*mu;
% aaa = mu11^2*eep;
% eeu = mean(diag(Dm2/LEN1));
% eeccc = 1-(1-mu11*eeu)^2;
% stdeeem(mm) = aaa/eeccc;
end
err = err_hampel1;
th_err = stdeeem;