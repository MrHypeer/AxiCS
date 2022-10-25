t = fft2c(I);
figure(2), imshow(abs(t),[]);
[mm,nn] = size(t);
ii = 0;
for kk = 1:1:mm
    for jj = 1:1:nn
        if abs(t(kk,jj))<20
            t(kk,jj) = 0;
            ii = ii+1;
        end 
    end
end 
tt = ifft2c(t);
ii = (mm*nn - ii)/mm/nn
subplot(3,2,1);
imshow(log(1+abs(t)),[]);
subplot(3,2,2);
imshow(abs(tt),[]);
 subplot(3,2,3), imshow(log(1+abs(fft2c(I))),[]);
 subplot(3,2,4), imshow(abs(I),[]);
  subplot(3,2,5), imshow(log(1+abs(fft2c(U))),[]);
 subplot(3,2,6), imshow(abs(U),[]);
 