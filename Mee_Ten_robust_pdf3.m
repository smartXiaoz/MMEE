close all
clear all
clc
 
N = 200;
C = 10
nn = 1 : N;
xx = C * (nn - N) / N;
y = zeros(N);
 
 
Nsam = 2000;
samp = randn(2, Nsam);
 
for iii = 1 : N
    iii;
    for jjj = 1 : N
        x(1) = C * (iii - N / 2) / N;
        x(2) = C * (jjj - N / 2) / N;
        
        for kkk = 1 : Nsam
            y(iii, jjj) = y(iii, jjj) + exp(-0.5 * (norm(x.' - samp(:, kkk)) ^2));
        end
    end
end
 
figure; mesh(xx, xx, y)
figure; imagesc(xx, xx, y)
