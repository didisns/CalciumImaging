function [ks,ksstd, area, areastd, cdf] = kstestGauss (flag, inData, lowLimit, highLimit, nFreedom) 
%
% This function computes the probability that the input data set (inData) is extracted from a
% normal distribution by using a Kolmogorov-Smirnov test against a standard Gaussian CDF
% lowLimit and highLimit are the -inf and +inf limits of the entire data
% set expressed in SD units

% GMR April 6, 2018.

% KS based on the kstest2.m function for KS testing

% d: K-S statistique D corresponding to the maximum difference between the sample distribution 
% and the equivalent normal distribution
% p: probability that the inData data set is extracted from a normal distribution.
% flag: used for initialisation of the function. To be used in the calling
% program to initialize the unitary Gauss distribution and others
% parameters

persistent binN binN0
persistent gaussCDF 
persistent binLimits x2
persistent n j

newLen = length(inData);
% everything in the following if is executed only at the first call of the
% function 
if flag,
    % First, compute the cumulative distribution of the reference gaussian
    % Mean and SD are both set = 1. The fluorescence data for each pixel will
    % be normalized in the same way and therefore we will be able to use
    % the same Gaussian CDF for all pixels. Pretty smart, I'd say.
    
    % First, let's generate the x axis that will be used for the Gaussian and the pixel CDFs
    % The following takes 0.06 sec for 2000 frames.
    x2 = (lowLimit:0.01:highLimit)/sqrt(2);
    binLimits = (lowLimit:0.01:highLimit);
    gaussCDF = 1/2*erfc(-x2);
    gaussCDF = gaussCDF';
    binN = length(binLimits);
    binN0 = floor(binN/2);
    j = (1:101)';         % used for the Q function

    % n1, 2 are the degree of freedom of the data set. The normal CDF has
    % the same value as the data. If nFreedom is fixed there is no need to compute
    % lambda at each iteration
    if nFreedom > 0
        n1     =  nFreedom;
        n2     =  nFreedom;
        n      =  n1 * n2 /(n1 + n2);
    end
end
% Standardization of the sampled variable (data) with rejection of NaN values
inDataClean = inData(isfinite(inData));
if nFreedom == 0
    nFreedom = length(inDataClean);      % used for the computation of the KS statistics
    n  =  nFreedom * nFreedom / (2 * nFreedom);
end
meanDI = mean(inDataClean);
stdDI = std(inDataClean);
normData = (inDataClean - meanDI) / std(inDataClean);

if ~isempty(inDataClean)
    % Continue only if the array is not empty
    % Calculate the empirical (i.e., sample) CDF.
    binCounts  =  histc (normData , binLimits, 1);
    sumCounts  =  cumsum(binCounts)./sum(binCounts);
    cdf = sumCounts;

    % Compute the metric of the comparison.
    % I have impletend 2x2 criteria:
    % 1) KS comparison 
    % 07 04 2018 Added a new criteria: just looking at the difference of
    % the CDFs for x larger than zero (right tail: Ca transient!)

    deltaCDF  =  abs(sumCounts - gaussCDF);
    
    % beging the computation of the metrics
    KSstatistic   =  max(deltaCDF);
    % 2-sided test (default).
    %  Use the asymptotic Q-function to approximate the 2-sided P-value.
    lambda =  max((sqrt(n) + 0.12 + 0.11/sqrt(n)) * KSstatistic, 0);
    p  =  2 * sum((-1).^(j-1).*exp(-2*lambda*lambda*j.^2));
    p  =  min(max(p, 0), 1);
    ks = (1 - p);     % 1-p gaussian distribution
    ksstd = ks*stdDI ;
%       d = max(abs(sumCounts(binN0:binN) - gaussCDF(binN0:binN)));
       %d = sum(abs(sumCounts(binN0:binN) - gaussCDF(binN0:binN)));
    area = sum(abs(sumCounts - gaussCDF));
    areastd = stdDI * area;
else
    % set p to one if there are no valid data in the set
    ks = 0;
    ksstd = 0;
    area = 0;
    areastd = 0;
    cdf(1:binN) = 0;
end
