function [t] = time(d2,d3,p,n,E)
	s4 = n - 1 - d2 - d3;
	if (d2 + d3 == n - 1)
		t = 0;
	else
		recVal = 0;
		for i = 0:min(d3-d2,s4)
			k = 0;
			if i == 0
				k = 1;
			end
			m = min(1 + d2, s4 - i);
			sumVal = 0;
			for j = k:m
				sumVal = sumVal + (pd(i,j,p,d2,d3,s4) * E((d2+1)+i,(d3+1)+j));
			end
			recVal = recVal + sumVal;
		end
		t = (1 / (1 - pd(0, 0, p, d2, d3, s4))) * (1 + recVal);
	end
end

