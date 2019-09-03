% crop a image%
 cd 'G:\NUCL 655 Project\10\frame1-50'
 clc
 clear all
 close all

pic = imread('1001.jpg');
imshow(pic);

[x,y] = ginput(2); %??????????ginput????????????
pic_1 = imcrop(pic,[x(1),y(1),abs(x(1)-x(2)),abs(y(1)-y(2))]);

imname = dir('.\*.jpg');
im_num = length(imname);    %

for i=1:im_num
    I = imread(imname(i).name,'jpg');
    pic_1 = imcrop(I,[x(1),y(1),abs(x(1)-x(2)),abs(y(1)-y(2))]);
    filename=[num2str(i),'.jpg'];
    imwrite(pic_1,filename,'jpg');
end
