function parameter = TestParameterMetrix(m, t, s, bp)
    raw = MonteCarlo(m, t, s);
    data = raw(:, s+1:t+s+1, :);
    [cbp, css, y, z, ut] = ChowTest(data(:,:,1)', bp, 4, 1, 199);
    parameter = GenerateParameterSeries(y, z, ut, 199);
end
