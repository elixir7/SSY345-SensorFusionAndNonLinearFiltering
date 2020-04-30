function [x, y] = samplesToState(measurements, s1, s2)
    t1 = measurements(1,:);
    t2 = measurements(2,:);
    
    x = (tan(t1)*s1(1) - tan(t2)*s2(1) + s2(2) - s1(2)) ./ (tan(t1) - tan(t2));
    y = tan(t1).*(x - s1(1)) + s1(2);
end