function [val,pos] = Regularize(product,Kin)
%   Reference:Needell D£¬Vershynin R. Uniform uncertainty principle and
%   signal recovery via regularized orthogonal matching pursuit. 
%   Foundations of Computational Mathematics, 2009,9(3): 317-334.  
    productabs = abs(product);
    [productdes,indexproductdes] = sort(productabs,'descend');
    for ii = length(productdes):-1:1
        if productdes(ii)>1e-6
            break;
        end
    end
    %Identify:Choose a set J of the K biggest coordinates
    if ii>=Kin
        J = indexproductdes(1:Kin);
        Jval = productdes(1:Kin);
        K = Kin;
    else%or all of its nonzero coordinates,whichever is smaller
        J = indexproductdes(1:ii);
        Jval = productdes(1:ii);
        K = ii;
    end
    %Regularize:Among all subsets J0¡ÊJ with comparable coordinates
    MaxE = -1;
    for kk = 1:K
        J0_tmp = zeros(1,K);iJ0 = 1;
        J0_tmp(iJ0) = J(kk);
        Energy = Jval(kk)^2;
        for mm = kk+1:K
            if Jval(kk)<2*Jval(mm)
                iJ0 = iJ0 + 1;
                J0_tmp(iJ0) = J(mm);
                Energy = Energy + Jval(mm)^2;
            else
                break;
            end
        end
        if Energy>MaxE
            J0 = J0_tmp(1:iJ0);
            MaxE = Energy;
        end
    end
    if MaxE<1e-10
        pos=0;
    else
    pos =J0;
    end
    val = productabs(J0);
end