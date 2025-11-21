function [Out] = Tucker_ARE(X,X1,known)
    % Missing data
    Nway = size(X);
    knownSet = zeros(Nway);
    knownSet(known) = 1;
    Partialmiss = find(knownSet==0);
    
    Xk = zeros(Nway);
    Xk(known) = X(known);
    Xmiss = Xk;

    for idx = 1:4
        inputSet(:,:,(idx-1)*3+1:idx*3) = knownSet(:,:,:);
        inputData(:,:,(idx-1)*3+1:idx*3) = X(:,:,:);
    end
    KN = find(inputSet == 1);
    data = inputData(KN);
    Norder = size(inputSet);
    
    % TMac Parameter
    opts.alpha_adj = 0;
    opts.rank_min = [1 1 1 1];
    opts.rank_max = [100,100,12];
    
    maxrank=[256 256 12]; 
    init_rank=[1 1 1];
    
    Addrank=addrank_Tucker;    
    mid_ssim=zeros(1,3);
    ssim_mid=cell(1,3);
    ssim_out=zeros(1,3);
    key=0;
    for k=1:4
        for  i = 1:3
            ranktube = init_rank;
            Xinter=interpolation(Xmiss,2);
            Out.ssim_inter{k,i}=ssim(Xinter,X1);
            for j=1:length(Addrank{i})   % add offset
                ranktube(i) = init_rank(i)+Addrank{i}(j);
                if init_rank(i)==maxrank(i) && k~=1
                    key=1;
                    continue
                end
                ranktube(i) = min(ranktube(i),maxrank(i));
                if ranktube(i)>0
                  for m=1:3
                    [X_inc,Y_inc,Out_inc] = TMac_ar(data,KN,Norder,ranktube,opts);
                    Mrec = zeros(Norder);
                    for mrec = 1:3
                        Mrec = Mrec+Out_inc.alpha(mrec)*Fold(X_inc{mrec}*Y_inc{mrec},Norder,mrec);
                    end
                    ssim_out(m)=ssim(Mrec(:,:,1:3),Xinter);
                  end
                  ssim_mid{i}(j)=mean(ssim_out(:));
                else
                  ssim_mid{i}(j)=0;
                end
            end
            if key == 0
                [mid_ssim(i), index] = max(ssim_mid{i}(:));
                init_rank(i) = min((init_rank(i)+Addrank{i}(index)),maxrank(i));
                disp(['epoch_',num2str(k), '_RankItem_',num2str(i),'_ssim_',num2str(mid_ssim(i)),'_Rank_', num2str(init_rank)])
            else
                switch i
                    case 1
                       disp(['Rank has reached the maximum value','_Rank_ ',num2str(init_rank),'_ssim_',num2str(mid_ssim(N-1))])
                    otherwise
                      disp(['Rank has reached the maximum value','_Rank_ ',num2str(init_rank),'_ssim_',num2str(mid_ssim(i-1))])
                end
            end
            key=0;
            [X_inc1,Y_inc1,Out_inc1] = TMac_ar(data,KN,Norder,init_rank,opts);
            Mrec1 = zeros(Norder);
            for mrec1 = 1:3
                Mrec1 = Mrec1+Out_inc1.alpha(mrec1)*Fold(X_inc1{mrec1}*Y_inc1{mrec1},Norder,mrec1);
            end    

            Out.rank{k,i}=init_rank;

            [X_inc2,Y_inc2,Out_inc2] = TMac_ar(data,KN,Norder,init_rank,opts);
            Mrec2 = zeros(Norder);
            for mrec2 = 1:3
                Mrec2 = Mrec2+Out_inc2.alpha(mrec2)*Fold(X_inc2{mrec2}*Y_inc2{mrec2},Norder,mrec2);
            end

            T = zeros(Nway);
            for idx = 1:4
                midAssign(:,:,:)=Mrec2(:,:,(idx-1)*3+1:idx*3);
                T = T+midAssign;
            end
            T = T/4;
            T(known) = X(known);
            [OutT,OutSet]=RS4Tucker(X,T,Partialmiss,known,0.5);

            KN = find (OutSet==1);
            data = OutT(KN);
            Out.Set{k,i} = OutSet;
            Out.Tensor{k,i} = OutT;
        end
    end
       
    [X_inc3,Y_inc3,Out_inc3] = TMac_ar(data,KN,Norder,init_rank,opts);
    Mrec3 = zeros(Norder);
    for mrec3 = 1:3
        Mrec3 = Mrec3+Out_inc3.alpha(mrec3)*Fold(X_inc3{mrec3}*Y_inc3{mrec3},Norder,mrec3);
    end
    Z = zeros(Nway);
    for idx = 1:4
        temp(:,:,:)=Mrec3(:,:,(idx-1)*3+1:idx*3);
        Z = Z+temp;
    end
    Z = Z/4;
    Z(known) = X(known);
    Out.Xc = Z;
end
