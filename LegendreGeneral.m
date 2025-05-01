clear
N=10;
P=5;

%linearly transforming x data
%xVec = 2*pi*rand(N,1);
xVec=linspace(0,2*pi,N);
yVec = sin(xVec)';
xVec= (1/pi)*xVec-1; % linear transofrmation from [0,2pi] to [-1,1]

%for loop that runs the size of the data
%A=[];
xData=linspace(-1,1,200);

for n=1:N
     for m=1:P+1                           
         A(n,m)=legendreP(m-1,xVec(n));

     end
end
A
betaVec= inv(A'*A)*A'*yVec;

for m= 1:length(xData)
    x=xData(m);
    sum=0;
    for n=1:P+1
        sum=sum+betaVec(n)*legendreP(n-1,x);
    end
    yData(m)=sum;
end

clf;
plot(xVec,yVec,'.','MarkerSize',15,'Color',[0.2 0.6 0.7]); hold on

plot(xData,yData,'-','LineWidth',2,'Color',[0.3 0.3 0.9]); hold on

