load('CUBIC_Students_Linear_Least_Squares.mat');

xVec = DATA(:,1);
yVec = DATA(:,2);
xVecSquared =zeros(157,1);
xVecCubed =zeros(157,1);
%for loop that runs the size of the data
for n=1:size(xVec)
   % squares the data 
    s=DATA(n,1)*DATA(n,1);
    xVecSquared(n,1)=s;
    %cubes the data
    c=s*DATA(n,1);
    xVecCubed(n,1)=c;
end
oneVector = ones(size(xVec));
A=[oneVector xVec xVecSquared xVecCubed];

betaVec= inv(A'*A)*A'*yVec;

plot(xVec,yVec,'.','MarkerSize',15,'Color',[0.2 0.6 0.7]); hold on

xData=linspace(-2,4,100);
%vectors of equally spaced points xData=-2:0.1:4;
yData=betaVec(1)+betaVec(2)*xData+betaVec(3)*xData.^2+betaVec(4)*xData.^3;
plot(xData,yData,'-','LineWidth',4,'Color',[0.3 0.3 0.9])