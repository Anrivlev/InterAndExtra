function[A,W]=interpolation(X,Y,C)
% A - function value (restored)
% X - matrix of arguments, X(i,:) - point i
% Y - column of values
% C - vector of arguments (value of function is not known)
for j=1:size(X,2)
    minj=min(min(X(:,j)),C(j));
    maxj=max(max(X(:,j)),C(j));
    X(:,j)=(X(:,j)-minj)/(maxj-minj);
    C(j)=(C(j)-minj)/(maxj-minj);
end
% X is now dimensionless ( 0<=X(i,j)<=1 for all i's and j's)
minY=min(Y);
maxY=max(Y);
Y=(Y-minY)./(maxY-minY);
% Y is now dimensionless too
% W is matrix of metric uncertainty
W=zeros(size(X,1),size(X,1));
for i=1:size(X,1)
    for j=1:size(X,1)
        if i==j
            W(i,j)=sum((X(i,:)-C).^2);
        else
            W(i,j)=sum((X(i,:)-C(1,:)).*(X(j,:)-C(1,:)));
        end
    end
end
% W is filled with values
Flag=0;
for i=1:size(X,1)
    if X(i,:)==C(1,:)
        Flag=i;
        break;
    end
end
if Flag==0
    A = dot(W\(ones(size(X,1),1)),Y)/dot(W\(ones(size(X,1),1)),(ones(size(X,1),1)));
else
    A=Y(Flag,1);
end
A = A*(maxY-minY)+minY;