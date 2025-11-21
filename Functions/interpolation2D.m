function [Y]=interpolation2D(X,n)
    I1 = 2;J1 = 2;                 % KA parameters
    R=256;C=256;
    Xmiss=X;
     X1=zeros(size(Xmiss));
     
     switch n
         case 2
          for channel = 1:1
          IMin0 = Xmiss(:,:,channel);
          IMin0= double(IMin0);
          Mask{channel} = double(~(IMin0==0));
          X1(:,:,channel) = interp_linear2D(IMin0,~Mask{channel});
          end
         otherwise
          X1=X1;
     end
     Y=X1;
end