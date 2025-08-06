function phase = phi(n,m,k,params,excludeOrder)
	phase = angle(eqn.c(n,k,params,excludeOrder) * ...
		conj(eqn.c(m,k,params,excludeOrder)));
end