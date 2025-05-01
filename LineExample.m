%Loads Data <-- 100x2 Matrix (1st col x, 2nd col y)
load('LINE_Students_Linear_Least_Squares.mat')

% in this example we want to find the best fit line: 
% y=b0+b1x


%vector of x
xVec =DATA(:,1);

%vector of y 
yVec= DATA(:,2);

oneVector = ones(size(xVec));

A=[oneVector xVec];

%Beta Coefficients for the model function
betaVec= inv(A'*A)*A'*yVec
