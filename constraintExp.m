function value = constraintExp(xIn,xLim)

x = max(0,xIn-xLim);

value = exp(x) - 1 - x;

end