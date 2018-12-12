function [ outputStr ] = extractBetween(text,leftStr,rightStr)
% Function extractBetween 
% This function extract a substring from text contained between the leftStr
% and rightStr strings.
% If the delimiter strings are not present in text the fuction returns an empty string

k1 = strfind (text, leftStr);
k2 = strfind (text, rightStr);

if isempty(k1),
    % the left string is not contained in text. Return with empty output
    outputStr = [];
    return
end

% find the first occurrency of leftStr
firstK1 = k1(1);
% find the first occurency of rightStr that follows firstK1
numRight = length(k2);
i = 1;
while (k2(i)<=firstK1) && i<=numRight
    % increment i
    i = i + 1;
end
if (k2(i) > firstK1),
    firstK2 = k2(i);
else 
    % if k2 firstK2 is not identified return empty string
    outputStr = [];
    return
end

outputStr = text (firstK1:firstK2);

