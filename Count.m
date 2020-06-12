function [CountOfCbp, CountOfCss] = Count(result)
    CountOfCbp = 0;
    CountOfCss = 0;
    for i = result'
        if i(1) > 0.5
            CountOfCbp = CountOfCbp + 1;
        end
        if i(2) > 0.5
            CountOfCss = CountOfCss + 1;
        end
    end
end
