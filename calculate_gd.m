function [d, Gd] = calculate_gd(N,Q)
    d_g(:,1) = zeros(Q-2,1);
    Gd_ = N*(N-1)/2;
    
    n = 1:N-1;
    
    % Calculate Uniques
    for d = 1:(n(end)/2)
        for i = 1:n(end-1-2*(d-1))
            d_g(:,end+1) = D(i, i+d, d, (3:Q)-1);
            d_g(:,end+1) = - d_g(:,end);
            Gd_(end+1) = 1;
            Gd_(end+1) = 1;
        end
    end
    
    % Calculate Repeaters
    for d = 1:ceil(n(end)/2-1)
        for j = 1:n(end-2*d)
            for i = 1:n(end-2*d-(j-1))
                d_g(:,end+1) = D(i, i+d, d+j, (3:Q)-1);
                d_g(:,end+1) = - d_g(:,end);
                Gd_(end+1) = 2;
                Gd_(end+1) = 2;
            end
        end
        
    end

    d = d_g; 
    Gd = Gd_;
end

function d = D(n,m,k,q)
    d = -(n+k).^q + (n).^q + (m+k).^q - m.^q;
end