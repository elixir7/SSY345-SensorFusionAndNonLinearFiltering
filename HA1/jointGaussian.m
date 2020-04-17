function [mu, Sigma] = jointGaussian(mu_x, sigma2_x, sigma2_r)
    %jointGaussian calculates the joint Gaussian density as defined
    %in problem 1.3a. 
    %
    %Input
    %   MU_X        Expected value of x
    %   SIGMA2_X    Covariance of x
    %   SIGMA2_R    Covariance of the noise r
    %
    %Output
    %   MU          Mean of joint density 
    %   SIGMA       Covariance of joint density


    %Your code here
    [mu, Sigma] = affineGaussianTransform([mu_x; 0], [sigma2_x 0; 0 sigma2_r], [1 0; 1 1], 0);

end