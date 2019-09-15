clc,clear all,close all
filename=('test_file.jpg')
img=imread(filename);
figure(1)
imshow(img)
img2=im2bw(img,150/225);
figure(2)
imshow(img2)
