function Mvv = mhampel(en,LEN)
L = LEN;

u = 0.6745 * en / median(abs(en - median(en)));
aa = 2;
bb = 4;
cc = 8;
w5 = zeros(L,1);
jw5 = zeros(1,L);
for jj = 1 : L
    if abs(u(jj)) <= aa
        w5(jj) = 1;
        Mvv(1,jj) = en(jj);        
    elseif (aa < abs(u(jj))) && (abs(u(jj)) <=bb)
        w5(jj) = aa / abs(u(jj));
        Mvv(1,jj) = en(jj) * w5(jj);       
    elseif (bb < abs(u(jj))) && (abs(u(jj)) <=cc)
        w5(jj) = (aa / abs(u(jj))) * ((cc-abs(u(jj)))/(cc-bb));
         Mvv(1,jj) = en(jj) * w5(jj);      
    elseif abs(u(jj)) > cc
        w5(jj) = 0;
        Mvv(1,jj) = 0;
    end
end
