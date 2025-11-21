function [X0] = iniMatrics(N,Nway2)
    IL =1;
    for k = 1:N-1
        dimL(k) = IL*Nway2(k);
        dimR(k) = prod(Nway2)/dimL(k);
        %
        X0{k} = randn(dimL(k),dimR(k));
        %
        IL = dimL(k);
        X0{k} = bsxfun(@rdivide,X0{k},sqrt(sum(X0{k}.^2,1)));

    end   
end