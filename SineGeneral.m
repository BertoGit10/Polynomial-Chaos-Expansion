N=50;

xVec = 2*pi*rand(500,1);
yVec = xVec.^3;
%for loop that runs the size of the data
A=[];
xData=linspace(0,2*pi,100);

for n=1:N
   A=[A sin((n*xVec)/2)];
    
end

betaVec= inv(A'*A)*A'*yVec;

for m= 1:length(xData)
    x=xData(m);
    sum=0;
    for n=1:N
        sum=sum+betaVec(n)*sin((x*n)/2);
    end
    yData(m)=sum;
end

clf;
plot(xVec,yVec,'.','MarkerSize',15,'Color',[0.2 0.6 0.7]); hold on

plot(xData,yData,'-','LineWidth',2,'Color',[0.3 0.3 0.9])