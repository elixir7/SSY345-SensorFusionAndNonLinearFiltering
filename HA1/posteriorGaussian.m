function [mu, sigma2] = posteriorGaussian(mu_x, sigma2_x, y, sigma2_r)
    %posteriorGaussian performs a single scalar measurement update with a
    %measurement model which is simply "y = x + noise".
    %
    %Input
    %   MU_P            The mean of the (Gaussian) prior density.
    %   SIGMA2_P        The variance of the (Gaussian) prior density.
    %   SIGMA2_R        The variance of the measurement noise.
    %   Y               The given measurement.
    %
    %Output
    %   MU              The mean of the (Gaussian) posterior distribution
    %   SIGMA2          The variance of the (Gaussian) posterior distribution

    % Calculate variance by expanding the exponentials and grouping 
    % terms of x together.
    % Then match these against our expected dist., which is a 
    % normal dist. And we get the following:
    sigma2 = sigma2_x * sigma2_r / (sigma2_x + sigma2_r);
    mu = sigma2_x/(sigma2_x + sigma2_r) * y + sigma2_r/(sigma2_x + sigma2_r) * mu_x;

end