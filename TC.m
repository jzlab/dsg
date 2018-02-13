function [ means ] = TC(angle, alpha,beta,gamma,phi)
%Gives outputs of Von Mises Dist

means = alpha + beta.*exp(gamma.*(cos(angle - phi) - 1));

end

