function deltaD = deltaD_mat(N,Q)
    arguments
        N (1,1) double {mustBePositive}
        Q (1,1) double {mustBePositive}
    end
	deltaD = zeros(N-1, N-1, N-1, Q);
	for q = 1:Q
		for k = 1:N-1
			for n=1:N-k
				for m=1:N-k

					deltaD(n,m,k,q) = eqn.deltaD(n,m,k,q);

				end
			end
		end
	end
end