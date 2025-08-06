function A = Amplitude(n,m,k,excludeOrder, params)
	A = 4*abs(ComputeResidualComplexWeight(n,k,params,excludeOrder) * ...
		conj(ComputeResidualComplexWeight(m,k,params,excludeOrder)));
end