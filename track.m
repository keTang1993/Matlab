clc
clear all
close all

% input a value   %标定

    prompt = 'Enter pix value:';
    opts.Interpreter = 'tex';
    title = 'Input';
    dims = [1 35];
    definput = {''};
    answer = inputdlg(prompt,title,dims,definput,opts);
    clear title
%} 
pix=str2num(answer{1});         % delta x

cd 'G:\NUCL 655 Project\3\video'   %进入视频文件夹
imname = dir('.\*.jpg');            %识别所有JPG格式的照片
im_num = length(imname);            %图片的总数
p=imread('background.jpg');               % background image
p1=imread('1001.jpg');                    % provessing image
fg1=rgb2gray(p);                 % 彩色到灰白to gray picture
fg2=rgb2gray(p1);                % 
D=abs(fg2-fg1);                  % subtract the background image
figure,imshow(D)
t=30;                           % boundary criteria，灰度值    
D(D<=t)=0;                       % 小于阈值的时候设为纯白
D(D>t)=255;                      % 大于阈值的时候设为纯黑

D= bwareaopen(D,100);            % remove objects with pixels less than 100
D= xor(bwareaopen(D,100),  bwareaopen(D,500));       % remove objects smaller than 30 pixels and larger than 500 pixels.
figure,imshow(D);

d=strel('disk',20);              % build structure element, to remove the inside white area of bubble, need to change base on bubble size
BW=imclose(D,d);                 % smooth the boundary, get rid of the rag,and fill the blank in the bubble
BW=~BW;                          % bubble black, background white
figure, imshow(~(BW-BW))
B=bwboundaries(BW);
hold on

%loop start here
h = waitbar(0);                 % 进度条
video=VideoWriter('trajectory.avi');    
open(video);
for i=1:im_num-1                % to eliminate the background image
    I = imread(imname(i).name,'jpg');
    fg2=rgb2gray(I);            % 
    D=abs(fg2-fg1);
    t=30;                      % boundary criteria    
    D(D<=t)=0;                  %
    D(D>t)=255;
    D= bwareaopen(D,200);        % remove objects with pixels less than 50
    D= xor(bwareaopen(D,30),  bwareaopen(D,2000));
    d=strel('disk',5);         % build structure element, to remove the inside white area of bubble, need to change base on bubble size
    BW=imclose(D,d);            % smooth the boundary, get rid of the rag,and fill the blank in the bubble
    BW=~BW;                     % bubble black, background white
    B=bwboundaries(~BW);        % B is a cell martrix,represent the boundary coordinate,
    imshow(BW)
    Bo=B{1,1};
    [Cen]=bwlabel(~BW);                             % label connected components in 2D binary image
    
    x_axis=regionprops(Cen,'MajorAxisLength');      % Majoraxis       
    y_axis=regionprops(Cen,'MinorAxisLength');      % Minoraxis
    Centroid = regionprops(Cen,'Centroid');
    Area=regionprops(Cen,'Area');
    Orientation=regionprops(Cen,'Orientation');
    bubble_geo(i,1)=x_axis.MajorAxisLength;         % bubbble horizontal diameter
    bubble_geo(i,2)=y_axis.MinorAxisLength;         % bubble vertical diameter
    bubble_geo(i,3)=Area.Area;                      % bubble area
    bubble_geo(i,4)=Centroid.Centroid(1);           % bubble centroid x-axis position
    bubble_geo(i,5)=Centroid.Centroid(2);           % bubble centroid y-axis position
    bubble_geo(i,6)= (bubble_geo(i,1).^2.*(bubble_geo(i,2))).^(1/3);  % Equivalent diameter 
    bubble_geo(i,7)=Orientation.Orientation;
    
   %每隔一帧画一个质心点
    if mod(i,1)==0
        plot(bubble_geo(i,4),bubble_geo(i,5),'b.','LineWidth',2)        % plot centroid 
    else
    end
    print('-dpng','trajectory.png')
    I = imread('trajectory.png');
    writeVideo(video,I);
    
  %每隔30帧画一个气泡轮廓
%     if mod(i,30)==0
%        plot(Bo(:,2),Bo(:,1),'r','LineWidth',2)                          % plot contour
%     else
%     end

    
    waitbar(i/im_num,h)
end
close(video)


for t=1:i-1
vz(t)= (bubble_geo(t,5)- bubble_geo(t+1,5))*pix/0.01+100;   % vertical velocity in units mm/s, local downward flow velocity aboout 100mm/s
vx(t)= (bubble_geo(t+1,4)- bubble_geo(t,4))*pix/0.01;       %horizontal velocity in units mm/s, 
delta_x(t)=(bubble_geo(t,4)- bubble_geo(1,4))*pix;          % lateral migration diatance
delta_z(t)=(bubble_geo(1,5)- bubble_geo(t+1,5))*pix;        %vertical migration distance
end


%output vertical velocity %
r=linspace(0,0.63,63);    % from 0.01s to 1.13s

figure                       % creat a new figure
plot(r,vz);               % plot vetical velocity vs time
xticks(0.1:0.1:0.7)
xlabel('t [s]');
ylabel('Vertical velocity [mm/s]')
title('Bubble vertical velocity')
print('-dpng','vertical_velocity.png')

% output horizontal velocity
figure                       % creat a new figure
plot(r,vx);               % plot horizontal velocity vs time
xticks(0.1:0.1:0.7)
xlabel('t [s]');
ylabel('Horizontal velocity [mm/s]')
title('Bubble horizontal velocity')
print('-dpng','horizontal_velocity.png')

figure
plot(r,delta_x)
hold on
p=polyfit(r,delta_x,1);    %%%%%%%plot fitting curve%%%%%%%%%%
y=polyval(p,r);
plot(r,y)
if p(2)<0
    p(2)=abs(p(2));
text(.1*max(r),.95*max(delta_x),sprintf('$(x-x_i)=%.3f(t-t_i)-%.3f$',p), 'interpreter',  'latex', 'FontSize', 10)
else
    text(.1*max(r),.95*max(delta_x),sprintf('$(x-x_i)=%.3f(t-t_i)+%.3f$',p), 'interpreter',  'latex', 'FontSize', 10)
end

xticks(0.1:0.1:0.7)
xlabel('t [s]');
ylabel('Lateral migration distance [mm]')
legend('Displpacement','Fitting curve')
print('-dpng','Horizontal_velocity_fit.png')
%}
save test.mat
mean(vz)