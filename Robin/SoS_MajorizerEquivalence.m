clc; clear; close all;
%% Sequence Set Design parameters
N = 8;
L = 2;
M = N/L;
% rng(1);
rng('shuffle');

Q = 3;

theta_li = [];
m = (1:M).';
a_coeff = zeros(3,L);
for l = 1:L
    a_coeff(:,l) = 2*pi*1e3*rand(Q,1);
    theta_li(:,l) = a_coeff(1,l)+ a_coeff(2,l).*m + a_coeff(3,l).*m.^2;
end
theta_x = theta_li(:);
x = exp(1j*theta_x);

theta_y = 2*pi*1e3*rand(N,1);
xg_prime = exp(1j*theta_y);

l = 1;
a = a_coeff(:,l);
indcs = (1+(l-1)*M: l*M)';
xg = zeros(N,1);
xg(indcs) = x(indcs);

theta_s = zeros(N,1);
theta_s(indcs) = theta_x(indcs);

theta_sa = ones(N,1).*(-1j*log(1e-16));
theta_sa(indcs) = theta_x(indcs);
x_tmpa = exp(1j*theta_sa);

x_g = x;
x_g(indcs) = 0;

c = abs(xcorr(x,xg_prime)).^2;
s_c = sum(c);
sk = 0;
sk_t = 0;
sk_ta = 0;

epsln = 1e-16;
f_tst = 4;
switch f_tst
    case 1
        % \sum_{k}\bigg|\Tilde{\mathbf{x}}_{s,g}^{H}\boldsymbol{\alpha}_{k}+{\gamma}_{k}\bigg|^{2}
        for k = 1-N:N-1
            Uk = toeplitzC(N,k);
            alpha_k = Uk*xg_prime;
            gamma_k = x_g' * Uk * xg_prime;

            x_tmp = zeros(N,1);
            x_t1 = exp(1j*theta_s);
            x_tmp(indcs) = x_t1(indcs);
            sk_t = sk_t + abs(x_tmp'*alpha_k + gamma_k)^2;
            sk   = sk   + abs(xg'   *alpha_k + gamma_k)^2;   % Error-free
            %             sk = sk + abs(exp(1j*2*pi*theta_s)'*alpha_k + gamma_k)^2; % Error
        end
        s_c
        sk
        sk_t

    case 2
        % \sum_{k}\bigg|\Tilde{\mathbf{x}}_{s,g}^{H}\boldsymbol{\alpha}_{k}+{\gamma}_{k}\bigg|^{2}
        for k = 1-N:N-1
            Uk = toeplitzC(N,k);
            alpha_k = Uk*xg_prime;
            gamma_k = x_g' * Uk * xg_prime;

            sk_ta = sk_ta + abs(x_tmpa'*alpha_k + gamma_k)^2;
            sk   = sk   + abs(xg'   *alpha_k + gamma_k)^2;   % Error-free
            %             sk = sk + abs(exp(1j*2*pi*theta_s)'*alpha_k + gamma_k)^2; % Error
        end
        s_c
        sk
        sk_ta

    case 3
        % \sum_{k}\bigg\{\bigg[\sum_{n}\psi_{nk}\bigg(\cos(\theta_{n}-\beta_{nk})\bigg)+{\gamma_{kr}}\bigg]^{2}
        sk = 0;
        for k = 1-N:N-1
            Uk = toeplitzC(N,k);
            tmp = Uk*xg_prime;
            alpha_k = zeros(N,1);
            alpha_k(indcs) = tmp(indcs);
            gamma_k = x_g' * Uk * xg_prime;

            z = 0;

            for n = 1:N
                theta_nk = theta_s(n);
                z = z + exp(-1i*theta_nk)*alpha_k(n);
            end

            sk = sk + abs(z + gamma_k)^2;
        end
        s_c
        sk

    case 4
        % \sum_{k}\bigg\{\bigg[\sum_{n}\psi_{nk}\bigg(\cos(\theta_{n}-\beta_{nk})\bigg)+{\gamma_{kr}}\bigg]^{2}
        for k = 1-N:N-1
            Uk = toeplitzC(N,k);
%             alpha_k = Uk*xg_prime;
            tmp = Uk*xg_prime;
            alpha_k = zeros(N,1);
            alpha_k(indcs) = tmp(indcs);
            gamma_k = x_g' * Uk * xg_prime;

            zr = 0;
            zi = 0;
            for n = 1:N
                psi_nk = sqrt(real(alpha_k(n)).^2 + imag(alpha_k(n)).^2) + epsln;
                theta_nk = theta_s(n);

                cos_beta_nk = (real(alpha_k(n))/psi_nk);
                sin_beta_nk = (imag(alpha_k(n))/psi_nk);
                rt = cos_beta_nk*cos(theta_nk) + sin_beta_nk*sin(theta_nk);
                it = sin_beta_nk*cos(theta_nk) - cos_beta_nk*sin(theta_nk);

                zr = zr + rt;
                zi = zi + it;
            end
            zri = zr + real(gamma_k);
            zki = zi + imag(gamma_k);

            znk = (zri^2 + zki^2);
            sk = sk + znk;
        end
        s_c
        sk

    otherwise
        % do nothing
end

%% -----------------------------------------  %
% ---------- Support functions -------------- %
%  -----------------------------------------  %

% Toeplitz function
function [Uk] = toeplitzC(N, k)
Uk = zeros(N);
for i = 1:N
    for j = 1:N
        if j-i == k
            Uk(i,j) = 1;
        else
            Uk(i,j) = 0;
        end
    end
end
end