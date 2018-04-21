function g = sigmoid(z)
%SIGMOID Compute sigmoid function
%   g = SIGMOID(z) computes the sigmoid of z.

% You need to return the following variables correctly 
dim1 = size(z, 1);
dim2 = size(z, 2);

g = zeros(dim1, dim2);

% ====================== YOUR CODE HERE ======================
% Instructions: Compute the sigmoid of each value of z (z can be a matrix,
%               vector or scalar).

for i=1:dim1
    for j=1:dim2
      g(i, j) = 1/(1+exp(-z(i, j)));
    endfor
endfor


% =============================================================

end
