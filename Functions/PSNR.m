function PSNR=PSNR(X,Y)
 %%%%%%%%%
  if size(X,3)~=1
      test=rgb2ycbcr(X);
      org=rgb2ycbcr(Y);
      
      Y1=test(:,:,1);
      Y2=org(:,:,1);
      Y1=double(Y1);
      Y2=double(Y2);
  else
      Y1=double(X);
      Y2=double(Y);
  end
  if nargin<2
      D=Y1;
  else
      if any(size(Y1)~=size(Y2))
          error('size is not equal');
      end
      D=Y1-Y2;
  end
  
  MSE=sum(D(:).*D(:))/numel(Y2);
  PSNR=10*log10(255^2/MSE);
end