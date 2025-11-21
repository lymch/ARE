function Out =tcu_left_contraction(G, C, order, n)

    N =size (n,2); 
    Grihgt =G {n(1)}; 
    for i =1:N-1
        Grihgt =tensor_contraction (Grihgt, G{n(i+1)}, 2*i+2, 1);    
    end
    C =permute (C,order);
    
    Grihgt =tensor_contraction (Grihgt, C, [3:2:5], [4:5]);  
    m = [2:3];
    n = [1,4:7];
    tempG =permute (Grihgt,[m,n]); 
    Nway = size(tempG);
    Out =reshape (tempG,prod(Nway(1:2)),prod(Nway(3:7))); 
end
