function V = logsigmoid(Q, theta, sigma, Qmax)
%LOGSIGMOID Summary of this function goes here
%   Detailed explanation goes here
    V = theta - sigma .* log(Qmax ./ Q - 1);
end

