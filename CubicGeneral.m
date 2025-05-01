load('CUBIC_Students_Linear_Least_Squares.mat');
N=8;

xVec = DATA(:,1);
yVec = DATA(:,2);
%for loop that runs the size of the data
oneVector = ones(size(xVec));
A=[oneVector];
xData=linspace(-2,4,100);



for n=1:N
   A=[A xVec.^n];
    
end

betaVec= inv(A'*A)*A'*yVec;

for m= 1:length(xData)
    x=xData(m);
    sum=0;
    for n=0:N
        sum=sum+betaVec(n+1)*x.^n;
    end
    yData(m)=sum;
end

clf;
plot(xVec,yVec,'.','MarkerSize',15,'Color',[0.2 0.6 0.7]); hold on

plot(xData,yData,'-','LineWidth',2,'Color',[0.3 0.3 0.9])