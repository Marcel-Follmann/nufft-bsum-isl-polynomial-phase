function value = deltaD(n, m, k, order)
	% Introduced in Eqn. (7)
	value = eqn.d(n,k,order) - eqn.d(m,k,order);
end