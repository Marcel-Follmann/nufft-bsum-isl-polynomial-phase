function amplitude = A(n,m,k,params, excludeOrder)
	amplitude = 4*abs(eqn.c(n,k,params,excludeOrder) * ...
		conj(eqn.c(m,k,params,excludeOrder)));
end