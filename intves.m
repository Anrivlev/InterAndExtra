function[A,W]=intves(X,Y,C)
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
A = dot(W\(ones(size(X,1),1)),Y)/dot(W\(ones(size(X,1),1)),(ones(size(X,1),1)));
wk=zeros(1,size(X,2));
for l=1:size(X,2)
    Xnew=X;
    Xnew(:,l)=[];
    Cnew=C;
    Cnew(l)=[];
    Wnew=zeros(size(Xnew,1),size(Xnew,1));
  for i=1:size(Xnew,1)
      for j=1:size(Xnew,1)
          Wnew(i,j)=sum((Xnew(i,:)-Cnew(1,:)).*(Xnew(j,:)-Cnew(1,:)));
      end
  end
  Wnew=pinv(Wnew);
  Ynew=dot(Wnew*(ones(size(Xnew,1),1)),Y)/dot(Wnew*(ones(size(Xnew,1),1)),(ones(size(Xnew,1),1)));
  wk(1,l)=(Ynew-A)^2;
end
w=zeros(1,size(X,2));
for i=1:size(X,2)
w(1,i)=(size(X,2)*wk(1,i))/(sum(wk(1,:)));
end


W=zeros(size(X,1),size(X,1));
for i=1:size(X,1)
    for j=1:size(X,1)
        if i==j
            W(i,j)=sum(w(1,:).*(X(i,:)-C).^2);
        else
            W(i,j)=sum(w(1,:).*(X(i,:)-C(1,:)).*(X(j,:)-C(1,:)));
        end
    end
end
A = dot(W\(ones(size(X,1),1)),Y)/dot(W\(ones(size(X,1),1)),(ones(size(X,1),1)));

A = A*(maxY-minY)+minY;