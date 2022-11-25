function p = InSight(z)
% InSight takes a vector input that contains 
% the positions of two particles. Checks if
% the line between them intersects the line
% segment x = -200, -400 <= y <= 0.
grad = (z(5)-z(2))/(z(4)-z(1));
y = grad.*(-200 - z(1)) + z(2);
if -400 <= y && y <= 0 
    p = false;
else 
    p = true;
end