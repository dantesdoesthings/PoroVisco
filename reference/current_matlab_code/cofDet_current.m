function [ RSquare adjRSquare ] = cofDet(expData, fitData, varargin)
%cofDet Calculates Rsquare coefficient of determination, and ajusted
%rsquare statistic for expdata and fit. Based on info in wiki article.
%   Varargin should be the num of regressors not including the constant!

n = length(expData);

if length(fitData) ~= n
    disp('*** lengths do not match ***')
    return
else
end

SSE = sum((expData-fitData).^2);
SST = sum((expData-mean(expData)).^2);
%SSR = sum((fitData-mean(expData)).^2);
RSquare = 1 - SSE/SST;

if nargin == 3
    p = varargin{1};
    adjRSquare = 1-(1-RSquare)*(n-1)/(n-p-1);
end

end

