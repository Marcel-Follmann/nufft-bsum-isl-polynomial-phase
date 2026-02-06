M = 1e5;                            % number of NU source points
x = 2*pi*rand(M,1);                 % points in a 2pi-periodic domain
c = randn(M,1)+1i*randn(M,1);       % iid random complex data (row or col vec)
N = 2e5;                            % how many desired Fourier modes?
opts = finufft_options();
opts.debug = 2;
epsilon = 1e-12;

w = ceil(log10(1/epsilon))+1;
beta = 2.30*w;
sigma = 2;

[q, p, r] = integer_factorization(sigma, N, w);
n = 2^q*3^p*5^r;
alpha = pi*w/n;
h = 2*pi/w;

f = finufft1d2(x,+1,1e-12,c, opts);     % do it (takes around 0.02 sec)

function [q,p,r] = integer_factorization(sigma, N, w)
    integer = max(sigma*N, 2*w);
    qmax = ceil(log(integer)/log(2));
    pmax = ceil(log(integer)/log(3));
    rmax = ceil(log(integer)/log(5));

    best = inf;
    best_exp = [0 0 0];
    
    for q = 0:qmax
        v2 = 2^q;
        if v2 >= best, break; end
    
        for p = 0:pmax
            v23 = v2 * 3^p;
            if v23 >= best, break; end
    
            for r = 0:rmax
                v = v23 * 5^r;
                if v >= integer && v < best
                    best = v;
                    best_exp = [q p r];
                end
                if v >= best, break; end
            end
        end
    end

    q = best_exp(1);
    p = best_exp(2);
    r = best_exp(3);
end