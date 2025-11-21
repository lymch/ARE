function [Y]=interpolation(X,n)
     Xmiss=X;
     X1=zeros(size(Xmiss));
     
     switch n
         case 1
          for channel = 1:3
          IMin0 = Xmiss(:,:,channel);
          IMin0= double(IMin0);
          Mask{channel} = double(~(IMin0==0));
          X1(:,:,channel) = interp_nearest(IMin0,~Mask{channel});
          end
         case 2
          for channel = 1:3
          IMin0 = Xmiss(:,:,channel);
          IMin0= double(IMin0);
          Mask{channel} = double(~(IMin0==0));
          X1(:,:,channel) = interp_linear(IMin0,~Mask{channel});
          end
         case 3
          for channel = 1:3
          IMin0 = Xmiss(:,:,channel);
          IMin0= double(IMin0);
          Mask{channel} = double(~(IMin0==0));
          X1(:,:,channel) = interp_d(IMin0,~Mask{channel});
          end
         case 4
          for channel = 1:3
          IMin0 = Xmiss(:,:,channel);
          IMin0= double(IMin0);
          Mask{channel} = double(~(IMin0==0));
          X1(:,:,channel) = interp_cubic(IMin0,~Mask{channel});
          end
         otherwise
          X1=X1;
     end
     Y=X1;
end