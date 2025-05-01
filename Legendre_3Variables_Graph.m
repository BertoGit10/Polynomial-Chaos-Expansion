clear
P=7; % amount of polynomials we are using 

%linearly transforming x data
xVec= legendre_root(P+1); % xVec will be the P polynomial roots

%---------------------------------------------------
%creates our alpha matrix which has all the combinations indexes
ct=0;
for a1=0:P
    for a2=0:P
        for a3=0:P
            if a1+a2+a3<=P
                ct=ct+1;
                alphaMATRIX(ct,1)=a1;
                alphaMATRIX(ct,2)=a2;
                alphaMATRIX(ct,3)=a3;
            end
        end
    end
end
%alphaMATRIX
%----------------------------------------------------
%create root matrix, with all our legendre polynomial roots
ct=0;
for r1=1:length(xVec)
    for r2=1:length(xVec)
        for r3=1:length(xVec)
            ct=ct+1;
            rootMATRIX(ct,1)=xVec(r1);
            rootMATRIX(ct,2)=xVec(r2);
            rootMATRIX(ct,3)=xVec(r3);
        end
    end
end

%---------------------------------------------------
polynomialTest=(1/8)*((3*(0.5*rootMATRIX(:,1)+0.5).^2+1).*(3*(0.5*rootMATRIX(:,2)+0.5).^2+1).*(3*(0.5*rootMATRIX(:,3)+0.5).^2+1));
ishigami= sin(pi.*rootMATRIX(:,1))    +(7*sin((pi.*rootMATRIX(:,2))).^2)+    (0.1*((pi.*rootMATRIX(:,3)).^4).*sin(pi.*rootMATRIX(:,1)));
yVec= ishigami;
%-----------------------------------------------
%Matrix A sets up the overconstrained linear system
for n=1:length(rootMATRIX(:,1))
     for m=1:length(alphaMATRIX(:,1))                          
         A(n,m)=legendreP(alphaMATRIX(m,1),rootMATRIX(n,1))*legendreP(alphaMATRIX(m,2),rootMATRIX(n,2))*legendreP(alphaMATRIX(m,3),rootMATRIX(n,3));
     end
end

%--------------------------------------------------
betaVec = (A' * A) \ (A' * yVec);
%------------------------------------------------------------------

totalVariance=0;

for n=2:length(betaVec)
    test=betaVec(n)*betaVec(n);
    first= 1/(2*alphaMATRIX(n,1)+1);
    second= 1/(2*alphaMATRIX(n,2)+1);
    third= 1/(2*alphaMATRIX(n,3)+1);
    totalVariance= totalVariance + test*first*second*third;
end 
totalVariance;





%-----------------------------------------------------------------%
%STORING S1 rows that matter ( second and third index are 0)
S1=1;
for rowNumber=2:length(alphaMATRIX)
    if alphaMATRIX(rowNumber,2)+alphaMATRIX(rowNumber,3)==0
       indexVectorS1(S1)=rowNumber;
       S1=S1+1;
       end 
end
%indexVectorS1         


%calculating partial variance of S1
sumS1=0;
clacS1=0;
for ct=1:length(indexVectorS1)
    calcS1=betaVec(indexVectorS1(ct))^2*(1/(2*alphaMATRIX(indexVectorS1(ct),1)+1));
    sumS1 = sumS1 +calcS1;
end
sumS1= sumS1/totalVariance;
sumS1 %partial variance of S1
%-----------------------------------------------------------------%


S2=1;
for rowNumber=2:length(alphaMATRIX)
    if alphaMATRIX(rowNumber,1)+alphaMATRIX(rowNumber,3)==0
       indexVectorS2(S2)=rowNumber;
       S2=S2+1;
       end 
end
%indexVectorS2       


%calculating partial variance of S2
sumS2=0;
clacS2=0;
for ct=1:length(indexVectorS2)
calcS2=betaVec(indexVectorS2(ct))^2*(1/(2*alphaMATRIX(indexVectorS2(ct),2)+1));
sumS2 = sumS2 +calcS2;
end
sumS2= sumS2/totalVariance;
sumS2 %partial variance of S2
%-----------------------------------------------------------------%

S3=1;
for rowNumber=2:length(alphaMATRIX)
    if alphaMATRIX(rowNumber,1)+alphaMATRIX(rowNumber,2)==0
       indexVectorS3(S3)=rowNumber;
       S3=S3+1;
       end 
end
%indexVectorS3       


%calculating partial variance of S3
sumS3=0;
clacS3=0;
for ct=1:length(indexVectorS3)
    calcS3=betaVec(indexVectorS3(ct))^2*(1/(2*alphaMATRIX(indexVectorS3(ct),3)+1));
    sumS3 = sumS3 +calcS3;
end
sumS3 =sumS3/totalVariance;

sumS3 %partial variance of S3
%---------------------------------------------------------------------
%second order index for parameters 1 and 2
% we are looing at alpha mat rows where the first and second parameters are
% not zero, and the 3rd parameter is
S12=1;
for rowNumber=2:length(alphaMATRIX)
    if alphaMATRIX(rowNumber,1)~=0 && alphaMATRIX(rowNumber,2)~=0 && alphaMATRIX(rowNumber,3)==0
       indexVectorS12(S12)=rowNumber;
       S12=S12+1;
    end 
end
sumS12=0;
clacS12=0;

for ct=1:length(indexVectorS12)
    calcS12=betaVec(indexVectorS12(ct))^2*(1/(2*alphaMATRIX(indexVectorS12(ct),1)+1))*(1/(2*alphaMATRIX(indexVectorS12(ct),2)+1));
    sumS12 = sumS12 +calcS12;
end

sumS12 =sumS12/totalVariance;

sumS12

%-----------------------------------------------------------------------
%second order index of parameters 1 and 3
S13=1;
for rowNumber=2:length(alphaMATRIX)
    if alphaMATRIX(rowNumber,1)~=0 && alphaMATRIX(rowNumber,2)==0 && alphaMATRIX(rowNumber,3)~=0
       indexVectorS13(S13)=rowNumber;
       S13=S13+1;
    end 
end
sumS13=0;
clacS13=0;

for ct=1:length(indexVectorS13)
    calcS13=betaVec(indexVectorS13(ct))^2*(1/(2*alphaMATRIX(indexVectorS13(ct),1)+1))*(1/(2*alphaMATRIX(indexVectorS13(ct),3)+1));
    sumS13 = sumS13 +calcS13;
end

sumS13 =sumS13/totalVariance;

sumS13

%-----------------------------------------------------------------------

%second order index of parameters 2 and 3
S23=1;
for rowNumber=2:length(alphaMATRIX)
    if alphaMATRIX(rowNumber,1)==0 && alphaMATRIX(rowNumber,2)~=0 && alphaMATRIX(rowNumber,3)~=0
       indexVectorS23(S23)=rowNumber;
       S23=S23+1;
    end 
end
sumS23=0;
clacS23=0;

for ct=1:length(indexVectorS23)
    calcS23=betaVec(indexVectorS23(ct))^2*(1/(2*alphaMATRIX(indexVectorS23(ct),2)+1))*(1/(2*alphaMATRIX(indexVectorS23(ct),3)+1));
    sumS23 = sumS23 +calcS23;
end

sumS23 =sumS23/totalVariance;

sumS23

%-----------------------------------------------------------------------
%third order index of parameters 1 and 2 and 3

S123=1;
for rowNumber=2:length(alphaMATRIX)
    if alphaMATRIX(rowNumber,1)~=0 && alphaMATRIX(rowNumber,2)~=0 && alphaMATRIX(rowNumber,3)~=0
       indexVectorS123(S123)=rowNumber;
       S123=S123+1;
    end 
end
sumS123=0;
clacS123=0;

for ct=1:length(indexVectorS123)
    calcS123=betaVec(indexVectorS123(ct))^2*(1/(2*alphaMATRIX(indexVectorS123(ct),1)+1)*(1/(2*alphaMATRIX(indexVectorS123(ct),2)+1))*(1/(2*alphaMATRIX(indexVectorS123(ct),3)+1)));
    sumS123 = sumS123 +calcS123;
end

sumS123 =sumS123/totalVariance;

sumS123

%------------------------------------------------------------------
% total for parameter 1


T1=1;
for rowNumber=2:length(alphaMATRIX)
    if alphaMATRIX(rowNumber,1)~=0 
       indexVectorT1(T1)=rowNumber;
       T1=T1+1;
    end 
end
sumT1=0;
clacT1=0;

for ct=1:length(indexVectorT1)
    calcT1=betaVec(indexVectorT1(ct))^2*(1/(2*alphaMATRIX(indexVectorT1(ct),1)+1)*(1/(2*alphaMATRIX(indexVectorT1(ct),2)+1))*(1/(2*alphaMATRIX(indexVectorT1(ct),3)+1)));
    sumT1 = sumT1 +calcT1;
end

sumT1/totalVariance


%----------------------------------------------------------------
% total for parameter 2


T2=1;
for rowNumber=2:length(alphaMATRIX)
    if alphaMATRIX(rowNumber,2)~=0 
       indexVectorT2(T2)=rowNumber;
       T2=T2+1;
    end 
end
sumT2=0;
clacT2=0;

for ct=1:length(indexVectorT2)
    calcT2=betaVec(indexVectorT2(ct))^2*(1/(2*alphaMATRIX(indexVectorT2(ct),1)+1)*(1/(2*alphaMATRIX(indexVectorT2(ct),2)+1))*(1/(2*alphaMATRIX(indexVectorT2(ct),3)+1)));
    sumT2 = sumT2 +calcT2;
end

sumT2/totalVariance

%----------------------------------------------------------------
% total for parameter 3


T3=1;
for rowNumber=2:length(alphaMATRIX)
    if alphaMATRIX(rowNumber,3)~=0 
       indexVectorT3(T3)=rowNumber;
       T3=T3+1;
    end 
end
sumT3=0;
clacT3=0;

for ct=1:length(indexVectorT3)
    calcT3=betaVec(indexVectorT3(ct))^2*(1/(2*alphaMATRIX(indexVectorT3(ct),1)+1)*(1/(2*alphaMATRIX(indexVectorT3(ct),2)+1))*(1/(2*alphaMATRIX(indexVectorT3(ct),3)+1)));
    sumT3 = sumT3 +calcT3;
end

sumT3/totalVariance


%----------------------------------------------------------------
% checking if it adds up to 1
% sumS1+sumS2+sumS3+sumS12+sumS13+sumS23+sumS123

% Visualization (Slice View)
% clf;
% figure(1)
% surf(X1,X2,gPC)
% xlabel('x1')
% ylabel('x2')
% zlabel('GPC')
% title('gPC Prediction')
% colorbar
% 
% figure(2)
% surf(X1,X2,abs(gPC-realData))
% xlabel('x1')
% ylabel('x2')
% zlabel('GPC')
% title('Error')
% colorbar
