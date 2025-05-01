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

yVec=rootMATRIX(:,1).^5 + rootMATRIX(:,2).^5;

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
x1Data=linspace(-1,1,10);
x2Data=linspace(-1,1,10);
[X1,X2]=meshgrid(x1Data,x2Data);

for m1= 1:length(x1Data)  %evaluating legendre polnomial series at test data points
    for m2 =1:length(x2Data)
    x1=x1Data(m1);
    x2=x2Data(m2);
    sum=0;
    for n=1:length(alphaMATRIX(:,1))
        sum=sum+betaVec(n)*legendreP(alphaMATRIX(n,1),x1)*legendreP(alphaMATRIX(n,2),x2);
    end
    gPC(m2,m1)=sum;%at points m1,m2 we are predicting a value which is given by sum
    realData(m2,m1)=(x1)^5+(x2)^5; %computes the real function we are looking for
    end
end
clf;
figure(1)
surf(X1,X2,gPC)
xlabel('x1')
ylabel('x2')
zlabel('gpC Prediction')
colorbar
figure(2)
surf(X1,X2,abs(gPC-realData))
xlabel('x1')
ylabel('x2')
zlabel('error')