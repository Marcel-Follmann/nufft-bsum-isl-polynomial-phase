close all; clearvars; clc;

N = 11;
Q = 5; 


[G, d_g] = eqn.G(zeros(Q,1), Q, eqn.deltaD_mat(N, Q), N, Q);
d_g(G == 0) = [];
G(G == 0) = [];
G(d_g<0) = [];
d_g(d_g<0) = [];


dn(:,1) = zeros(Q-3,1);
Gn = N*(N-1)/2;

n = 1:N-1;

% Calculate Uniques
for d = 1:(n(end)/2)
    for i = 1:n(end-1-2*(d-1))
        dn(:,end+1) = D(i, i+d, d, (3:Q-1));
        Gn(end+1) = 1;
    end
end

% Calculate Repeaters
for d = 1:ceil(n(end)/2-1)
    for j = 1:n(end-2*d)
        for i = 1:n(end-2*d-(j-1))
            dn(:,end+1) = D(i, i+d, d+j, (3:Q-1));
            Gn(end+1) = 2;
        end
    end
    
end

stem(d_g, G, 'LineWidth',2,'MarkerSize',8); hold on; 
plot(dn(2,:), 2*Gn, 'X', 'LineWidth',2, 'MarkerSize',15);

function d = D(n,m,k,q)
    d = -(n+k).^q + (n).^q + (m+k).^q - m.^q;
end