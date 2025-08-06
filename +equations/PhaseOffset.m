function phase = PhaseOffset(n,m,k,params,excludeOrder)
	phase = angle(ComputeResidualComplexWeight(n,k,params,excludeOrder) * ...
		conj(ComputeResidualComplexWeight(m,k,params,excludeOrder)));
end