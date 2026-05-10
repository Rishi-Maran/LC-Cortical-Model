function Q = sigmoid(V, theta, sigma, Qmax)
%SIGMOID Summary of this function goes here
%   Detailed explanation goes here
    Q = Qmax ./ (1 + exp(-(V - theta)./sigma));
end

