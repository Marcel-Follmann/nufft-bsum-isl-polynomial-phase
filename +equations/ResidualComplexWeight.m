function residual = ResidualComplexWeight(n,k, params, excludeOrder)
	% exp(j sum_{q=1 and q \neq p)^Q a_q d_{n,k}^{(q)}
	nOrder = numel(params);
	
	phaseSum = 0;
	for iOrder = 1:nOrder
		if iOrder == excludeOrder, continue; end

		phaseSum = phaseSum + ...
			params(iOrder) * ComputeShiftDifference(n,k,iOrder);
	end

	residual = exp(1i*phaseSum);
end