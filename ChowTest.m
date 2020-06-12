function [Chi_Cbp, Chi_Css, Y12, Z12, UT12] = ChowTest(Variables, BreakDate, LagOrder, InterceptIndicator, BtsRepNums)
    [t, k] = size(Variables);
    TargetData = lagmatrix(Variables, 1: LagOrder);

    Y1 = Variables(LagOrder+1: (BreakDate-1), :)';
    Z1 = TargetData(LagOrder+1: (BreakDate-1), :)';
    T1 = length(Y1);
    if InterceptIndicator == 1
        Z1 = [ones(1,T1); Z1];
    end

    Y2 = Variables(LagOrder+BreakDate: end, :)';
    Z2 = TargetData(LagOrder+BreakDate: t, :)';
    T2 = length(Z2);
    if InterceptIndicator == 1
        Z2 = [ones(1,T2); Z2];
    end


    Y12 = Variables(LagOrder+1: t, :)';
    Z12 = TargetData(LagOrder+1: end, :)';
    T12 = length(Z12);
    if InterceptIndicator == 1
        Z12 = [ones(1,T12); Z12];
    end

    %@ is the syntax of anoymous function
    LSEstimator = @(Y, Z) Y*Z'*inv(Z*Z');
    GetUT = @(Y, Z)  Y - LSEstimator(Y, Z) * Z;

    UT1 = GetUT(Y1, Z1);
    UT2 = GetUT(Y2, Z2);
    UT12 = GetUT(Y12, Z12);


    SigTilde1 = 1 / T1 * UT1 * UT1';
    SigTilde2 = 1 / T2 * UT2 * UT2';


    temp1 = UT12(:, 1: T1);
    temp2 = UT12(:, (T12-T2+1): T12);
    SigTilde12 = 1 ./ (T1+T2) * (temp1 * temp1' + temp2 * temp2');

    loglM1 = -(T1/2) * log(det(SigTilde1)) - (T2/2) * log(det(SigTilde2));
    loglM2 = -(T1+T2) / 2 * log(det(SigTilde12));


    [t, K] = size(Variables);
    % OutPut Valut Of ChowTest
    ChowBreakPoint = 2*(loglM1 - loglM2);
    DegreeOfFreedom = K^2 * LagOrder + K + (K * (K+1))/2;
    Chi_Cbp = 1 - chi2cdf(ChowBreakPoint, DegreeOfFreedom);
    ChowSampleSplit = (T1 + T2) * (log(det(SigTilde12)) - log(det((T1+T2)^-1 * (T1*SigTilde1 + T2*SigTilde2))));
    DegreeOfChowSampleSplit = K^2 * LagOrder + K;
    Chi_Css =  1 - chi2cdf(ChowSampleSplit, DegreeOfChowSampleSplit);
end
