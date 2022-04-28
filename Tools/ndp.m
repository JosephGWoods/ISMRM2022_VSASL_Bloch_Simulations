%% Function to return the number of non-zero decimal places
%
% num = ndp(x, maxNumDecimalPlaces)
%
% If x is an array, the function will find the largest number of decimal
% places.
%
% Written by Joseph G Woods, FMRIB, University of Oxford, 2018

function num = ndp(x, maxNumDecimalPlaces)

if ~exist('maxNumDecimalPlaces','var') || isempty(maxNumDecimalPlaces)
    maxNumDecimalPlaces = 12;
end

x = x(:);

% floating point errors are sometime created here, so only use the minimum
% maxNumDecimalPlaces that are necessary.
x = round( x - floor(x), maxNumDecimalPlaces);

if ~any(x~=0)
    num = 0;
else
    y = x*10.^(0:maxNumDecimalPlaces);
    
    num = zeros(size(y,1),1);
    for ii = 1:size(y,1)
        try    num(ii) = find(y(ii,:)==round(y(ii,:)),1) - 1;
        catch; num(ii) = 0; end
    end
    
    num = max(num);
end

end