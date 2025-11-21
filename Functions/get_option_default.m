function option=get_option_default(MissM,Mask,alg,opt)

option=[];

framenum=opt.framenum;

if ndims(MissM)==3&&framenum~=1
    option.nChannel=1;
else
    option.nChannel=3;
end

option.stopc=1e-5;
option.maxitr  = 500;
option.debug=0;

switch alg
    case 'PTRC'
        n_t=opt.n_t;
        N_t=length(n_t);
        L=ceil(N_t/2);
        p=sum(double(Mask(:)==1))/sum(double(Mask(:)));
        sk=[];
        for nn=1:N_t
            order=[nn:N_t 1:nn-1];
            M=reshape(MissM,prod(n_t(order(1:L))),[]);
            sk=[sk max(ceil((min(size(M)))*opt.eta*sqrt(p)),max(floor(sqrt(framenum)*2),3))];
        end
        option.d=ceil(N_t/2);
        option.beta=1/N_t*ones(N_t,1);
        
        option.lambda=1;
end
end