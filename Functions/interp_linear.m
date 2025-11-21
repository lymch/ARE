function im_r=interp_linear(im_n,c)
    [x,y]=find(c==0);
    [M,N]=size(im_n);
    [x1,y1]=meshgrid(1:M,1:N);
    temp = im_n(find(c==0));
    im_r=griddata(x,y,temp,x1,y1,'linear');  
    im_r=im_r';
    I=find(isnan(im_r)==1);
     im_r(I)=128;
end