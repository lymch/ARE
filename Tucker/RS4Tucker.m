function [Out,finalSet] = RS4Tucker(X,Assign,Ps,known,ratio)
    T = zeros(size(Assign));
    knownSet = T;
    knownSet(known) = 1;
    knownSet1 = knownSet;knownSet2 = knownSet;knownSet3 = knownSet;
    
    knownSet(known) = 1;
    finalSet = zeros([256 256 12]);
    finalSet(:,:,1:3) = knownSet(:,:,:);
    
    sumT = prod(size(T));
    num=size(Ps,1);
    numSelect = num*ratio;
    
    addOne=randsample(sumT,round(numSelect));
    addTwo=randsample(sumT,round(numSelect));
    addThree=randsample(sumT,round(numSelect));
    knownSet1(addOne) = 1;
    knownSet2(addTwo) = 1;
    knownSet3(addThree) = 1;

    finalSet(:,:,4:6) = knownSet1(:,:,:);
    finalSet(:,:,7:9) = knownSet2(:,:,:);
    finalSet(:,:,10:12) = knownSet3(:,:,:);

    Out(:,:,1:3) = X(:,:,:);
    for idx = 2:4
        Out(:,:,(idx-1)*3+1:idx*3) = Assign(:,:,:);
    end
   
end

