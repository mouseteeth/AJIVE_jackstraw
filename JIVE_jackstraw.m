function [F_feature,F_n,P] = JIVE_jackstraw(Y,x,k,nran, nsim);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Y       datablock m*n
%  x       CNS for test Joint component r_j*n; 
%          BSSindiv 
%  k       number of component can be 1:k
%  nran    number of random chosen row for each permutation
%  nsim    number of permutated satistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Y=D{1};x=outstruct.BSSindiv{1};
for i=1:length(x(:,1))
    x(i,:)=x(i,:)/norm(x(i,:));
end
n=length(Y(1,:));d=length(Y(:,1));
if length(nran)==0
nran=1;
end
if length(nsim)==0
nsim=1000;
end
X=[ones(n,1)/sqrt(n) x(k,:)'];
Y=Y';
B=X'*Y;
p=length(X(1,:));
y_hat=X*B;
y_bar=mean(Y);
SSE=sum((Y-y_hat).^2);
SSM=sum((y_hat- y_bar).^2);   
F_feature=(SSM/(p-1))./ (SSE/(n-p)); 
%%%% F_n %%%%%
F_n=[];
for i=1:(nsim/nran)
    Y_per=Y;
    index=randsample(d,nran);
   for j=1:length(index)
       Y_per(:,index(j))=Y_per(randperm(n),index(j));
   end
   B=X'*Y_per(:,index);
   y_hat=X*B;
   y_bar=mean(Y_per(:,index));
SSE=sum((Y_per(:,index)-y_hat).^2);
SSM=sum((y_hat- y_bar).^2);   
F_n=[F_n (SSM/(p-1))./ (SSE/(n-p))];    
end
P=[];
for i=1:d
P=[P mean(F_feature(i)<F_n)];
end









