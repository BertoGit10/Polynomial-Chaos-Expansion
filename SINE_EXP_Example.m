load('SINE_EXP_Linear_Least_Squares.mat')
%vector of x
xVec =DATA(:,1);

%vector of y 
yVec= DATA(:,2);

xVecSIN =zeros(157,1);
xVecEXP =zeros(157,1);
%for loop that runs the size of the data
for n=1:size(xVec)
   % Sines the data 
    s=DATA(n,1);
    xVecSIN(n,1)=sin(2*s);
    %exponeniates the data
    c=DATA(n,1);
    xVecEXP(n,1)=exp(c/4);
end
oneVector = ones(size(xVecSIN));
A=[oneVector xVecSIN xVecEXP];

betaVec= inv(A'*A)*A'*yVec;
plot(xVec,yVec,'.','MarkerSize',15,'Color',[0.2 0.6 0.7]); hold on

xData=linspace(0,20,100);
%vectors of equally spaced points xData=-2:0.1:4;
yData=betaVec(1)+betaVec(2)*sin(2*xData)+betaVec(3)*exp(xData/4);
plot(xData,yData,'-','LineWidth',4,'Color',[0.3 0.3 0.9])