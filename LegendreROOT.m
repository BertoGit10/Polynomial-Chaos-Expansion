clear
P=4; % amount of polynomials we are using 
N=P+1; % amount of enumerations 

%linearly transforming x data
xVec= legendre_root(P+1)'; % xVec will be the P polynomial roots
yVec = sin(exp(xVec));

%for loop that runs the size of the data
A=[];

for n=1:N
     for m=1:N                           
         A(n,m)=legendreP(m-1,xVec(n));

     end
end

betaVec= inv(A'*A)*A'*yVec;

xData=linspace(-1,1,200);
for m= 1:length(xData) %evaluating legendre polnomial series at test data points
    x=xData(m);
    sum=0;
    for n=1:P+1
        sum=sum+betaVec(n)*legendreP(n-1,x);
    end
    yData(m)=sum;
    realyData(m) = sin(exp(xData(m)));; %computes the real function we are looking for
end

clf;
plot(xVec,yVec','.','MarkerSize',15,'Color',[0.2 0.6 0.1]); hold on% input points using the roots

plot(xData,yData,'-','LineWidth',2,'Color',[0.3 0.3 0.9]); hold on % trained function
plot(xData,realyData,'-','LineWidth',2,'Color',[0,0,0]);

