function [Out] = SobelInitial(X,X1,known,ratio,Xk9D,it,iter)
 
%% Missing data
    Nway = size(X); N = numel(Nway);
    R=256;C=256;I1 = 2; J1 = 2;
    knownSet=zeros(Nway);
    knownSet(known)=1;
    knownSet3D = CastKet2Image(knownSet,R,C,I1,J1);
%     Xk = CastKet2Image(Xk9D,R,C,I1,J1); 

    known9D=known;
    known3D=find(knownSet3D~=0);
    init_known9D=known9D;
    Xkn = X(known9D);
    Partialmiss = find(knownSet3D==0);

    init_known=known3D;%输入的是9D的
    Xmiss = Xk9D;
    
    
    
    %% TMac-TT Algorithm        
    tol = 10^(-4);      % tol from Eq. (43)
    maxiter = 500;     % max iterations
    
    maxrank=[4 16 64 256 192 48 12 3]; %设定秩上界
    rankBound=[4 16 64 100 100 48 12 3]; %设定秩上界
    init_rank=ones(1,N-1);
       
    mid_ssim=zeros(1,N-1);
    ssim_mid=cell(1,N-1);
    ssim_out=zeros(1,3);

    n=1;
    key=0;
    Xinter=interpolation(Xmiss,2);
    
    pop_size = 50;
    population = round(sobol_initialization(pop_size, N-1, rankBound));
    for pop = 1:pop_size
        [Xpop,~,~] = TMacTT(Xkn,known9D,Nway,N,population(pop,:),tol,maxiter);
        P{pop} = CastKet2Image(Xpop,R,C,I1,J1);
        ref(pop)=ssim(P{pop},Xinter);
    end
    startRank = zeros(1,N-1);
    [sortedRef, sortIdx] = sort(ref, 'descend');
    startRank = population(sortIdx(1),:)
%     startRank = round(sum(population(sortIdx(1:10),:))/10)
    
    init_rank = startRank;

    %%小范围筛选
       for k=1:it
           Addrank=offSet(k); 
           for  i = 1:8     %1到8阶
                ssim_mid=cell(1,N-1);
                ranktube = init_rank;
                for j=1:7  %每一阶的秩+候选秩
                  ranktube(i) = init_rank(i)+Addrank(j);
                  ranktube(i) = min(ranktube(i),maxrank(i));
                  if ranktube(i)>0

                     [Xout,~,~] = TMacTT(Xkn,known9D,Nway,N,ranktube,tol,maxiter);
                      mid_Xout = CastKet2Image(Xout,R,C,I1,J1);
                      ssim_mid{i}(j)=ssim(mid_Xout,Xinter);
                      %disp(['第',num2str(k),'次', '第',num2str(i),'阶','第',num2str(j),'候秩', num2str(ranktube),'对应ssim值',num2str(ssim_mid{i}(j)),' 对应PSNR',num2str(psnr_mid{i}(j))])
                  else
                      ssim_mid{i}(j)=0;
                  end
                end
                [selectRelerr, index] = max(ssim_mid{i}(:));
                init_rank(i) = min((init_rank(i)+Addrank(index)),maxrank(i));
                disp(['第',num2str(k),'次', '第',num2str(i),'阶最大值','ssim',num2str(selectRelerr),' 秩', num2str(init_rank)])

                [Xmid,~,~] = TMacTT(Xkn,known9D,Nway,N,init_rank,tol,maxiter);
                Xmid1 = CastKet2Image(Xmid,R,C,I1,J1);       

                Out.rank{k,i}=init_rank;
                Out.XC{k,i}=Xmid1;
        %         Out.knownset{k,i}=knowN;
                Out.ssim{k,i}=ssim(Xmid1,X1);
                Out.ssiminter{k,i}=ssim(Xmid1,Xinter);
%                 Out.Xk{k,i}=Xk;
                Out.Xmiss=Xmiss;
        %         Out.diff{k,i}=a; 
                [Out.Xc,~,~] = TMacTT(Xkn,known9D,Nway,N,init_rank,tol,maxiter);
           end
       end

       
end

