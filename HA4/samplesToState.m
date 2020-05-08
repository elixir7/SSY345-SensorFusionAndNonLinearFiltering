function [x, y] = samplesToState(measurements, S)
    distance = measurements(1,:);
    angle = measurements(2,:);
    
    pos = S + distance .* [cos(angle); sin(angle)];
    x = pos(1,:);
    y = pos(2,:);
    
end