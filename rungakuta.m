function [eDelay, phiDelay, delayIndex] = rungakuta(a, h, t, M, phaseShift, delayIndex, eDelay, phiDelay)
    x = zeros(1, M);
	b = zeros(4, M);
    theta = rem(phaseShift + a(1) - phiDelay(delayIndex), 2.0 * pi);
	for i = 1 : 4
		for j = 0 : M
			if i == 0
				x(j) = a(j);
			if i == 1
				x(j) = a(j) + h * b(0)(j) / 2.0;
			if i == 2
				x(j) = a(j) + h * b(1)(j) / 2.0;
			if i == 3
				x(j) = a(j) + h * b(2)(j);
        	laser(x, b[i], theta);
		end
	end
	for i = 1: M
		a[i] += h * (b[0][i] + 2.0 * b[1][i] + 2.0 * b[2][i] + b[3][i]) / 6.0;
	end
    eDelay[delayIndex] = a[0];
    phiDelay[delayIndex] = a[1];
    delayIndex = (delayIndex + 1);
end
