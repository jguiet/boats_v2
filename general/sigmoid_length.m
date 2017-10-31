function sigmoid = sigmoid_length(x,x_50,slope)
%----------------------------------------------------------------------------------------
% sigmoid_length(x,x_50,slope)
% calculate values of sigmoid function with input values in terms of length
% sigmoidal relationship is assumed to be in terms of length
% convert to mass to length and calculate sigmoid values
%----------------------------------------------------------------------------------------
% x				mass array
% x_50			parameter value at which function value is 0.50
% slope			parameter value that determines rate of increase of function
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% assume that mass = delta1 * length^delta2
% convert to length
 delta1 = 0.01;
 delta2 = 3;
 y = (x/delta1).^(1/delta2);
 y_50 = (x_50/delta1).^(1/delta2);

% declare sigmoid values array as zeros
 sigmoid = zeros(1,length(y));
% calculate sigmoid values based on input parameters
 sigmoid = 1 ./ (1 + exp( - (slope * (y - y_50) )));
% bound between 0 and 1
 sigmoid = min(1,max(0,sigmoid));

%-----------------------------------------------------------------------------------------
% END OF FUNCTION