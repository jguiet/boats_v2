function sigmoid = sigmoid_mass(x,x_50,slope)
%----------------------------------------------------------------------------------------
% sigmoid_mass(x,x_50,slope)
% calculate values of sigmoid function with input values in terms of mass
%----------------------------------------------------------------------------------------
% x				mass array
% x_50			parameter value at which function value is 0.50
% slope			parameter value that determines rate of increase of function
%-----------------------------------------------------------------------------------------

% declare sigmoid values array as zeros 
 sigmoid = zeros(1,length(x));
% calculate sigmoid values based on input parameters
 sigmoid = 1 ./ (1 + exp( - (slope * (x - x_50) )));
% bound value between 0 and 1
 sigmoid = min(1,max(0,sigmoid));
 
%----------------------------------------------------------------------------------------
% END OF FUNCTION