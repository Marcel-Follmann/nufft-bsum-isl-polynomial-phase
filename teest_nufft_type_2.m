close all; clearvars; clc; 
% Barnett, A. H., Magland, J. F., & Klinteberg, L. af. (2018). A parallel non-uniform fast Fourier transform library based on an „exponential of semicircle" kernel (Version 2). arXiv. https://doi.org/10.48550/ARXIV.1808.06736

% Input params
N = 20;
k(:,1) = 0:N-1; 
f_k = zeros(size(k)); 
f_k(7) = 1;

% Output params
M = 1000;
x_j = sort(2*pi*rand(M, 1));





b_k = h*f_k./kernel_transform(k, alpha, beta);
b_l = ifftshift(ifft(b_k, n));

l(1,:) = 0:n-1;
tmp = periodic_kernel(l*h - x_j, alpha, beta);
c_j = sum(b_l.' .* tmp,2);

nexttile;
plot(k, f_k);
nexttile;
plot((x_j - pi) / (2*pi), real(c_j), '-o');



function x = trunctuate(x, N)
    M = ceil(numel(x)/2);

    x = x(M-N/2:M+N/2-1);

end

function phi = periodic_kernel(z, alpha, beta)
    x = mod(z, 2*pi);
    m(1,1,:) = [-1,0,1]*2*pi;
    phi = sum(kernel(x - m, alpha, beta),3);
end

function phi = kernel(z, alpha, beta)
    z = z/alpha;
    phi = zeros(size(z));
    idx = abs(z)<=1;
    phi(idx) = besseli(0, beta*sqrt(1-z(idx).^2))./besseli(0,beta);
end

function phi = kernel_transform(z, alpha, beta)
    z = z*alpha;
    tmp = sqrt(beta^2 - z.^2);
    phi = alpha*(2/besseli(0, beta)) * sinh(tmp)./tmp;
end

function n = calculate_dft_size(sigma, N, w)
    [q, p, r] = factorization(max(sigma*N, 2*w));
    n = 2^q*3^p*5^r;
end

function [q, p, r] = factorization(integer)
    factors = Inf;
    while(any(factors > 5))
        factors = factor(integer);
        integer = integer + 1;
    end
    p = sum(factors == 3);
    q = sum(factors == 2);
    r = sum(factors == 5);
end