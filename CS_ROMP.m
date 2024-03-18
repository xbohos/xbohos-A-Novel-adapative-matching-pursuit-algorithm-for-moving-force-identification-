function [ theta ] = CS_ROMP( y,A ,K)
    %% Parameters
    [M,N] = size(A);
    %% 1
    theta = zeros(N,1);
    At = zeros(M,3*K);
    Pos_theta = zeros(1,2*K);
    Index = 0;
    r_n = y;
    %Repeat the following steps K times(or until |I|>=2K)
    for ii=1:size(A,2)
        product = A'*r_n;
        [~,pos] = Regularize(product,K);
        if pos==0
            break
        end
        At(:,Index+1:Index+length(pos)) = A(:,pos);
        Pos_theta(Index+1:Index+length(pos)) = pos;
        if Index+length(pos)<=M
            Index = Index+length(pos);
        else
            break;
        end
        A(:,pos) = zeros(M,length(pos));
        theta_ls = inv(At(:,1:Index)'*At(:,1:Index))*At(:,1:Index)'*y;
        r_n = y - At(:,1:Index)*theta_ls;
        if norm(r_n)<1e-6%Repeat the steps until r=0
            break;
        end
        if Index>=2*K%or until |I|>=2K
            break;
        end
    end
    theta(Pos_theta(1:Index))=theta_ls;
end