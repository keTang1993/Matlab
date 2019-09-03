clc
clear all
close all
cd 'G:\NUCL 655 Project\19\video'
video_file='19.avi';
video=VideoReader(video_file);
frame_number=floor(video.Duration * video.FrameRate);
 h = waitbar(0,'Please wait...');
 
for i=1:frame_number
    image_name=strcat(num2str(i+1000));
    image_name=strcat(image_name,'.jpg');
    I=read(video,i);                              
    imwrite(I,image_name,'jpg');                  
    I=[];
waitbar(i/frame_number,h)
end
delete(h)