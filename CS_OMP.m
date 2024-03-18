function[x]=CS_OMP(y,A,t)
%11.22
    [y_rows,y_columns] = size(y);
    if y_rows<y_columns
        y=y';
    end
    [M,N]=size(A);
    x=zeros(N,1);
    t=round(t);
    At=zeros(M,t);
   B=zeros(1,t);
    r_n=y;
    for i=1:t
        neiji=A'*r_n;
        [~,pos] = max(abs(neiji));
        At(:,i)=A(:,pos);
        B(i)=pos;
        A(:,pos)=zeros(M,1);
        xt=inv(At(:,1:i)'*At(:,1:i))*At(:,1:i)'*y;
        r_n=y-At(:,1:i)*xt;
        if r_n<1e-6
            B=B(1:length(xt));
            break
        end  
    end
    x(B)=xt;