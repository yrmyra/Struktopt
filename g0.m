function value = g0(u,F,z)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    value=F'*u*(sum(z)/length(z))^(3/2);

end

