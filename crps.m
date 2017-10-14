function s = crps(x,x0)
% CRPS computes the continuously ranked probability score of the sample x
% with respect to the value x0

x = sort(x,'ascend'); % Sort the sample
n = length(x);
s = 0; % Initialise sum

if x(1)<x0 && x0<x(n) % If x0 is covered by the sample
    ind = find(x<=x0,1,'last');
    for i = 2:ind
        s = s + (x(i)-x(i-1))*((i-1)/n)^2;
    end
    s = (x0-x(ind))*ind/n + (x(ind+1)-x0)*(1 - (ind+1)/n);
    for i = ind+2:n
        s = s + (x(i)-x(i-1))*((1 - (i-1)/n))^2;
    end
elseif x0<=x(1) % If x0 is to the left of the sample
    s = s + (x(1)-x0);
    for i = 2:n
        s = s + (x(i)-x(i-1))*((1 - (i-1)/n))^2;
    end
elseif x0>=x(n) % If x0 is to the right of the sample
    for i = 2:n
        s = s + (x(i)-x(i-1))*((i-1)/n)^2;
    end
    s = s + (x0-x(n));
else
    return;
end

end