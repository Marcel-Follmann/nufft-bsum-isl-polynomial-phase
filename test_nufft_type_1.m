close all; clearvars; clc; 
% Barnett, A. H., Magland, J. F., & Klinteberg, L. af. (2018). A parallel non-uniform fast Fourier transform library based on an „exponential of semicircle" kernel (Version 2). arXiv. https://doi.org/10.48550/ARXIV.1808.06736

% Input params
x(:,1) = pi/2*[-1,-1/10, 1/10, 1];
c(:,1) = 1/2*ones(numel(x),1);
M = numel(x); 

% Output params
N = 200;
s(:,1) = 0:N-1;

sigma = 2;
epsilon = 1e-12;
w = ceil(log10(1/epsilon)) + 1;
beta = 2.3 * w;
n = calculate_dft_size(sigma, N, w);
alpha = pi*w/n;
h = 2*pi/n;


l(1,:) = 0:n-1;
b_l(:,1) = sum(c.*periodic_kernel(l*h - x, alpha, beta),1);
b_k = fftshift(fft(b_l));
k(:,1) = 0:n-1;

b_k = trunctuate(b_k, N);
k = trunctuate(k - (n)/2, N);

f_k = h*b_k./(kernel_transform(k, alpha, beta));



nexttile; 
plot(linspace(-1.2,1.2,100), kernel(linspace(-1.2,1.2,100)*alpha, alpha, 4));
title("Kernel")
nexttile; 
plot(linspace(0,n+ N/2,1000), log(abs(kernel_transform(linspace(0,n+N/2,1000), alpha, beta)))); hold on; 
yLim = ylim;
plot((n) * [1,1], yLim, 'k');
plot((n-N/2) * [1,1], yLim, 'k');
plot((N/2) * [1,1], yLim, 'k');
title("Kernel transform");
nexttile; 
stem(x, c);
title("Function");
nexttile;
plot(l*h - pi, b_l); 
title("Spreading")
nexttile; 
plot(k*h, b_k);
title("Transform");
nexttile;
plot(k*h, f_k, '-o');
title("Deconvolution");

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