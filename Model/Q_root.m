function Q=Q_root(nus,theta,sigma,Qmax)
    % Find solutions for the steady state equation
    Q = zeros(5, 1);
    Q(1) = fzero(@(Qe) Qe_root(nus,Qe, theta, sigma, Qmax), [5 6], optimset('TolFun',1e-50));
    Q(2) = Q(1);
    Q(4) = (1/nus(1,4))*(logsigmoid(Q(1),theta,sigma,Qmax) - ...
        (nus(1,1) + nus(1,2))*Q(1));
    Q(3) = sigmoid(nus(3,4)*Q(4) + nus(3,1)*Q(1),theta,sigma,Qmax);
    Q(5) = 1;
end