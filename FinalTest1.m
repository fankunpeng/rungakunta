function result = FinalTest1(m, t, s, bp)
    raw = MonteCarlo(m, t, s);
    data = raw(:, s+1:t+s+1, :);
    result = [];
    for i = 1:m
        [cbp, css, ut] = ChowTest(data(:,:,i)', bp, 4, 1, 199);
        result(end+1, :) = [cbp; css];
        residents(:, :, i) = ut;
    end
end
