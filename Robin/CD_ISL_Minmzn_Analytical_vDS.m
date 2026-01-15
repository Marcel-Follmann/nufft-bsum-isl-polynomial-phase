clc; clear; close all;

%% Parameters
N = 8;
L = 1;
Q = 3;
M = N / L;

maxIter = 10;
tol = 1e-9;

a_coeff = rand(Q, L);
[x, phi] = create_sequence(a_coeff, M, L);


%% Auto-correlation equivalence test
ISL = 0;
a2 = a_coeff(3);  % Extract initial a2
a2_i = a2;
r = zeros(N-1,1);
s1 = 0; s2 = 0;
for k = 1:N-1
    s1 = s1 + (N-k);
    for n = 1:N-k
        for np = 1:N-k
            if n < np
                psi_knnp = 2 * a2 * k * (n - np);
                s2 = s2 + 2*cos(psi_knnp);
            else
                continue
            end
        end
    end
end


ISL = s1+s2;
rm = xcorr(x);
ISLm = sum(abs(rm(N+1:2*N-1)).^2);

fprintf('ISL computed using my logic: %f \n', ISL);
fprintf('ISL computed using Matlab: %f \n', ISLm);



%% Support functions
function [x, phi_x] = create_sequence(a_coeff, M, L)
phi_i = [];
m = (1:M).';
for l = 1:L
    phi_i(:,l) = a_coeff(1,l) + a_coeff(2,l).*m + a_coeff(3,l).*m.^2;
end
phi_x = phi_i(:);
x = exp(1j*phi_x);
end