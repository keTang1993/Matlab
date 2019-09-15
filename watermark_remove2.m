clc,clear all, close all
I=imread('2.png');
figure(10), imshow(I);
t =I(590 : 620, 240 : 440, 1 : 3);%¹À¼ÆË®Ó¡ÇøÓò
figure(1), imshow(t);
 t1 = t(:,:,1);
 figure(2), imshow(t1);
 t2 = t(:, :, 2);
 figure(3), imshow(t2);
 t3 = t(:, :, 3);
 figure(4), imshow(t3);
 
 [m,n] = size(t2)
 for i = 1:m
    for j = 1:n
        if t2(i,j) >= 150
            t2(i,j) = 250;
        end
    end
 end
 
figure(5), imshow(t2);
 
for i = 1:m
    for j = 1:n
        if t3(i,j) >= 150
            t3(i,j) = 250;
        end
    end
end
 
figure(6), imshow(t3);
 
for i = 1:m
    for j = 1:n
        if t1(i,j) >= 150
            t1(i,j) = 255;
        end
    end
end
 
figure(7),imshow(t1);
 
for i = 1:m
    for j = 1:n
        t(i,j,1) = t1(i,j);
        t(i,j,2) = t2(i,j);
        t(i,j,3) = t3(i,j);
    end
end
figure(8), imshow(t);
 
for i = 1:31
    for j = 1:201
        I(i + 589, j + 239, 1:3)=t(i, j, 1:3);
    end
end
figure(9), imshow(I);