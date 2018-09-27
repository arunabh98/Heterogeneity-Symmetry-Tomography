function D2 = distfun(ZI,ZJ)
	D2 = zeros(size(ZJ, 1), 1);

	for i=1:size(ZJ, 1)
		D2(i) = abs(sum(ZI) - sum(ZJ(i, :)));
%         D2(i) = abs(sum(ZI) - sum(ZJ(i, :))).^2;
%         D2(i) = abs(sum((ZI - ZJ(i, :))^(1/3)));
	end
end
