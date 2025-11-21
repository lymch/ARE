function [Out] = TMac_SobelT_Patero(X,X1,known,SR)
 
%% Missing data
    Nway = size(X);
    R=256;C=256;I1 = 2; J1 = 2;
    Xk = zeros(Nway);
    Xk(known) = X(known);
    Xmiss = Xk;
    
    
    data = X(known);
    %% TMac-TT Algorithm
    opts.alpha_adj = 0;
    opts.rank_min = [1 1 1];
    opts.rank_max = [80,80,3];
    
    switch SR
        case 0.1
            rank_max = [25,25,3];
        case 0.2
            rank_max = [30,30,3];
        otherwise
            rank_max = [80,80,3];
    end  
    maxrank=[256 256 3]; %设定秩上界

    
    Xinter=interpolation(Xmiss,2);
    pop_size = 50;
    population = round(sobol_initialization(pop_size, 3, rank_max));
    for pop = 1:pop_size
        [Xpop,Ypop,Out_pop] = TMac_ar(data,known,Nway,population(pop,:),opts);
        Mrec = zeros(Nway);
        for mpop = 1:3
            Mpop = Mrec+Out_pop.alpha(mpop)*Fold(Xpop{mpop}*Ypop{mpop},Nway,mpop);
        end
        P{pop} = Mpop;
        ref(pop)=ssim(P{pop},Xinter);
    end
    [~, sortIdx] = sort(ref, 'descend');
    startRank = population(sortIdx(1),:);
    
    Out.startRank = startRank;
    init_rank = startRank;
    Addrank=addrank_Tucker;  
    
    mid_ssim=zeros(1,3);
    ssim_mid=cell(1,3);
    ssim_out=zeros(1,3);
%     key=0;
    
    %%自适应T分布变异
    initial_nu = 1;        % 初始自由度
    final_nu = 30;         % 最终自由度
    mutation_rate = 0.1;   % 变异率
    
    for k=1:4
        for  i = 1:3     %1到3阶
            nu = initial_nu + (final_nu - initial_nu) * (k / 4);
            ranktube = init_rank;
            Xinter=interpolation(Xmiss,2);
            Out.ssim_inter=ssim(Xinter,X1);
            for j=1:length(Addrank{i})   %每一阶的秩+候选秩
              ranktube(i) = init_rank(i)+Addrank{i}(j);
%               if init_rank(i)==maxrank(i) && k~=1
%                 key=1;
%                 continue
%               end
              ranktube(i) = min(ranktube(i),maxrank(i));
              if ranktube(i)>0
                  for m=1:3
                     [X_inc,Y_inc,Out_inc] = TMac_ar(data,known,Nway,ranktube,opts);
                     Mrec = zeros(Nway);
                     for mrec = 1:3
                         Mrec = Mrec+Out_inc.alpha(mrec)*Fold(X_inc{mrec}*Y_inc{mrec},Nway,mrec);
                     end
                     ssim_out(m)=ssim(Mrec,Xinter);
                  end
                  ssim_mid{i}(j)=mean(ssim_out(:));
              else
                  ssim_mid{i}(j)=0;
              end
            end

            [mid_ssim(i), index] = max(ssim_mid{i}(:));
            init_rank(i) = min((init_rank(i)+Addrank{i}(index)),maxrank(i));
            disp(['第',num2str(k),'次', '第',num2str(i),'阶最大值','ssim',num2str(mid_ssim(i)),' 秩', num2str(init_rank)])

            
            %变异
            disturbance = t_distribution_mutation(nu, mutation_rate);
            dis_rank = init_rank;
            dis_rank(i) = round(min((init_rank(i) + init_rank(i)*disturbance),maxrank(i)));
            if dis_rank(i)<0
                dis_rank(i) = 1;
            end
            
            [X_inc1,Y_inc1,Out_inc1] = TMac_ar(data,known,Nway,init_rank,opts);
            Mrec1 = zeros(Nway);
            for mrec1 = 1:3
                Mrec1 = Mrec1+Out_inc1.alpha(mrec1)*Fold(X_inc1{mrec1}*Y_inc1{mrec1},Nway,mrec1);
            end    
            opt_ssim = ssim(Mrec1,Xinter);
            
            [X_inc2,Y_inc2,Out_inc2] = TMac_ar(data,known,Nway,dis_rank,opts);
            Mrec2 = zeros(Nway);
            for mrec2 = 1:3
                Mrec2 = Mrec2+Out_inc2.alpha(mrec2)*Fold(X_inc2{mrec2}*Y_inc2{mrec2},Nway,mrec2);
            end    
            dis_ssim = ssim(Mrec2,Xinter);
            
            if opt_ssim >= dis_ssim
                t = 0;
                Out.ssiminter(k,i)=opt_ssim;
                Out.ssim(k,i)=ssim(Mrec1,X1);
                Out.XC{k,i}=Mrec1;
            else
                t = 1;
                init_rank(i) = dis_rank(i);
                Out.ssiminter(k,i)=dis_ssim;
                Out.ssim(k,i)=ssim(Mrec2,X1);
                Out.XC{k,i}=Mrec2;
            end
            Out.rank{k,i}=init_rank;
        end
    end
    
    [~,Idx] = sort(Out.ssiminter(:),'descend');
    for num = 1:3
        ParetoRank = Out.rank{Idx(num)};
        [X_Pareto,Y_Pareto,Out_Pareto] = TMac_ar(data,known,Nway,ParetoRank,opts);
        Mrec_Pareto = zeros(Nway);
        for i_Pareto = 1:3
            Mrec_Pareto = Mrec_Pareto+Out_Pareto.alpha(i_Pareto)*Fold(X_Pareto{i_Pareto}*Y_Pareto{i_Pareto},Nway,i_Pareto);
        end 
        ssimPareto(num) = ssim(Mrec_Pareto,Xinter);
    end

    [FinalX,FinalY,FinalOut] = TMac_ar(data,known,Nway,init_rank,opts);
    Mrec_Final = zeros(Nway);
    for i_Final = 1:3
        Mrec_Final = Mrec_Final+FinalOut.alpha(i_Final)*Fold(FinalX{i_Final}*FinalY{i_Final},Nway,i_Final);
    end 
    ssimPareto(num+1) = ssim(Mrec_Final,Xinter);
    
    [~,IdxPareto] = max(ssimPareto);
    if IdxPareto<= 5
        init_rank = Out.rank{Idx(IdxPareto)};
    end
    add = ssimPareto(IdxPareto);
    Out.ssiminter(k,i+1) = add;
    
    [~,cir] = sort(Out.ssiminter(:),'descend');
    if cir(1) == k*(i+1)
        [Xc,Yc,Outc] = TMac_ar(data,known,Nway,init_rank,opts);
        Mc = zeros(Nway);
        for ij = 1:3
            Mc = Mc+Outc.alpha(ij)*Fold(Xc{ij}*Yc{ij},Nway,ij);
        end
        Out.Xc = Mc;
    else
        Out.Xc = Out.XC{cir(1)};
    end

end
