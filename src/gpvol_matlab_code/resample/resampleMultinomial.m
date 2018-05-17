function [ indx ] = resampleMultinomial( w, varargin )

if isempty(varargin)
    M = length(w);
else
    M = varargin{1};
end
Q = cumsum(w);
Q(M)=1; % Just in case...

i=1;
while (i<=M),
    sampl = rand(1,1);  % (0,1]
    j=1;
    while (Q(j)<sampl),
        j=j+1;
    end;
    indx(i)=j;
    i=i+1;
end