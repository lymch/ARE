function Out =tcu_right_contraction(G, n)
    N =size (n,2); 
    Grihgt =G {n(1)}; 
    for i =1:N-1
        Grihgt =tensor_contraction (Grihgt, G{n(i+1)}, 2*i+2, 1);    
    end
    m = [2:2:6];
    n = [1:2:7,8];
    tempG =permute (Grihgt,[m,n]); 
    Nway = size(tempG);
    Out =reshape (tempG,prod(Nway(1:3)),prod(Nway(1:4))); 
end
