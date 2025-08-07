function value = A(n,m,k,params, excludeOrder)
	% Introduced in Eqn. (7)
	value = eqn.c(n,k,params,excludeOrder) * ...
		conj(eqn.c(m,k,params,excludeOrder));
end