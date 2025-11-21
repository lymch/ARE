function [ Out ] = SiLRTC_tt( x,D,T,N,th)
%SILRTT 此处显示有关此函数的摘要
%   此处显示详细说明
    thl = [0.01:0.005:0.025];       % Adjust this value for different thresholds
    for iter = 1:length(thl)
            x(D) = T(D); 
            x_unfold = cell(1,N-1);
            M  = cell(1,N-1);
            xk = cell(1,N-1);
            Nway2 = size(x);
            x_sum = zeros(Nway2);
            tol = 10^(-4);      % tol 
            k=1;
            relerr(1) = 1; 
            th = thl(iter);
            [~,ranktube{iter}] =  SVD_MPS_Rank_Estimation(T,th);
            beta = ranktube{iter};
            while relerr(k) > tol
                k=k+1;
                if k>500
                    break;
                end
                xlast = x;
                x_sum = x_sum * 0;
                x_unfold = iniMatrics(N,Nway2);

                for i = 1:N-1
                    x_unfold{i} = reshape(x,prod(Nway2(1:i)),prod(Nway2(i+1:N)));
                    [U,S,V] = svd_Truncate(x_unfold{i},ranktube{iter}(i),th);
                    M{i} = U*S*V;
                    xk{i} = reshape(M{i},Nway2);
                    x_sum = x_sum + (beta(i)*xk{i})/ sum(beta);
                end
                x = x_sum ;
                x(D) = T(D);
                rse = RSE(x(:),T(:));
                relerr(k) = abs(norm((x(:)-xlast(:)),'fro') / norm(x(:),'fro'));
                disp(['k=',num2str(k),'   rse:',num2str(rse), '   relerr:',num2str(relerr(k))])
            end
            Xc{iter} = x;
            RSe(iter) = rse;      
    end
    [Out.RSEmin, Idx] = min(RSe);
    Out.muf = thl(Idx);
    Out.MS = Xc{Idx};   

end

