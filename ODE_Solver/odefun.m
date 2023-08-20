function dy=odefun(r,y)
global m
dy=[y(2);m^2*y(1)-y(2)/r];
end