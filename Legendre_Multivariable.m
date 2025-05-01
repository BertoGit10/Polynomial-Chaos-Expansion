clear
P=3; % amount of polynomials we are using 

%linearly transforming x data
xVec= legendre_root(P+1); % xVec will be the P polynomial roots

%for loop that runs the size of the data

%---------------------------------------------------
%creates our alpha matrix which has all the combinations indexes
ct=0;
for a1=0:P
    for a2=0:P
        if a1+a2<=P
            ct=ct+1;
            alphaMATRIX(ct,1)=a1;
            alphaMATRIX(ct,2)=a2;
        end
    end
end

%alphaMATRIX


%----------------------------------------------------
%create root matrix, with all our legendre polynomial roots
ct=0;
for r1=1:length(xVec)
    for r2=1:length(xVec)
            ct=ct+1;
            rootMATRIX(ct,1)=xVec(r1);
            rootMATRIX(ct,2)=xVec(r2);
    end
end
%rootMATRIX
%---------------------------------------------------

yVec=rootMATRIX(:,1).^3 + rootMATRIX(:,2).^3;

%-----------------------------------------------
%Matrix A sets up the overconstrained linear system, which will help us
%find the coefficients for the legendre series
for n=1:length(rootMATRIX(:,1))
     for m=1:length(alphaMATRIX(:,1))                          
         A(n,m)=legendreP(alphaMATRIX(m,1),rootMATRIX(n,1))*legendreP(alphaMATRIX(m,2),rootMATRIX(n,2));

     end
end
%A
%--------------------------------------------------
betaVec= inv(A'*A)*A'*yVec;
betaVec;
% 
% %-----------------------------------------
% 
xData=linspace(-1,1,10);
ct=0;
for m1= 1:length(xData)  %evaluating legendre polnomial series at test data points
    for m2 =1:length(xData)
       
    x1=xData(m1);
    x2=xData(m2);
    sum=0;
    for n=1:length(alphaMATRIX(:,1))
        sum=sum+betaVec(n)*legendreP(alphaMATRIX(n,1),x1)*legendreP(alphaMATRIX(n,2),x2);
    end
    ct=ct+1;
    yData(ct)=sum;
    realyData(ct) =(x1)^3+(x2)^3; %computes the real function we are looking for
    end
end
realyData;
clf;
plot(realyData,'k-','LineWidth',4); hold on
plot(yData,'r-','LineWidth',2); % trained function

 
