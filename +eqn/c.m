function residual = c(n,k, params, excludeOrder)
	% exp(j sum_{q=1 and q \neq p)^Q a_q d_{n,k}^{(q)}
	nOrder = numel(params);
	
	phaseSum = 0;
	for iOrder = 1:nOrder-1
		if iOrder == excludeOrder, continue; end

		phaseSum = phaseSum + ...
			params(iOrder+1) * eqn.d(n,k,iOrder);
	end

	residual = exp(1i*phaseSum);
end