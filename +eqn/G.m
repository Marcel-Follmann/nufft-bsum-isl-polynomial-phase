function [g, d] = G(a, p, deltaD, N, Q)
	arguments
		a (:,1) double
		p (1,1) double {mustBePositive} % Selected parameter a_p
        deltaD (:,:,:,:) double
        N (1,1) double {mustBePositive}
        Q (1,1) double {mustBePositive}
    end
    
	d_p = deltaD(:,:,:,p);  
	dmin = min(d_p(:));  
	dmax = max(d_p(:));  
	d = dmin:dmax;
	id_offset = 1 - dmin;

	g = zeros(size(d));
	for k = 1:N-1
		for n = 1:N-k
			for m = 1:N-k

				id = deltaD(n,m,k,p);
				idx = id + id_offset;

				phase = 0;
				for q = 2:Q
					if q == p, continue; end
					phase = phase + a(q)*deltaD(n,m,k,q);
				end
				g(idx) = g(idx) + 2*exp(1j*phase);

			end
		end
    end

    % d(g == 0) = [];
    % g(g == 0) = [];
end