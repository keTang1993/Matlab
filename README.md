# Matlab

nothing important, I am just saving the code I used before here.

Bubble rising up trajectory:

![image](https://github.com/keTang1993/Matlab/blob/master/WeChat%20Image_20190818104016.png?raw=true)

代码，作者 李文强

clc
clear
vn_candidate=110;%%气泡上升过程处理个数限值
videolimit=10000;%%%飘浮气泡提取个数
Magnification=0.071315;%单位mm/pixel，
destinationini=['H:\Experiment\Measurement 20190121\'];
destinationpost=[destinationini,'postprocess code\'];%
path(path,destinationpost);%将程序文件添加到路径，这样不会进行报错
FileDirDiameter=dir('H:\Experiment\Measurement 20190121\Diameter 0.3-FPS 400-Main Stream-800 1280 plu*');%程序文件起始路径 
% filefoldname=cell(1,1);
filefoldname1=FileDirDiameter(1).name;
%% 气泡侧面面积限值设定
exp='(?<=Diameter)\s\d\.\d';
[s,e,strmatch]=regexp(filefoldname1,exp,'start','end','match');
if ~isempty(strmatch)
    Diameter=str2num(string(strmatch));
    if Diameter==0.3
       Arealimit=[300,2000];%直径为0.3mm时候气泡侧视面积限值
    end
    if Diameter==0.7
       Arealimit=[1000,6000];%直径为0.7mm时候气泡侧视面积限值
    end
end
if isempty(strmatch)
    exp='(?<=Diameter)\s\d';
    [s,e,strmatch]=regexp(filefoldname1,exp,'start','end','match');
    Diameter=str2num(string(strmatch));
    if Diameter==1
       Arealimit=[1000,10000];%直径为1mm时候气泡侧视面积限值
    end
     if Diameter==2
       Arealimit=[2000,10000];%直径为2mm时候气泡侧视面积限值
     end
    if Diameter==3
       Arealimit=[2000,10000];%直径为3mm时候气泡侧视面积限值
    end
    if Diameter==4
       Arealimit=[2000,10000];%直径为4mm时候气泡侧视面积限值
    end
     if Diameter==5
       Arealimit=[2000,20000];%直径为5mm时候气泡侧视面积限值
    end
end
%%
[s,e,strmatch]=regexp(filefoldname1,exp,'start','end','match');
if ~isempty(strmatch)
    Diameter=str2num(string(strmatch));
    if Diameter==0.3
       SurfacelifeArealimit=[40,700];%平均值76
          Surfaceliferect=[688,460,546,508]; 
    end
    if Diameter==0.7
       SurfacelifeArealimit=[100,1000];%平均值249
          Surfaceliferect=[665,307,600,550]; 
    end
end
if isempty(strmatch)
    exp='(?<=Diameter)\s\d';%判定dimater不带小数点
    [s,e,strmatch]=regexp(filefoldname1,exp,'start','end','match');
    Diameter=str2num(string(strmatch));
    if Diameter==1
           SurfacelifeArealimit=[300,3000];%平均值337,416,443
           if ~isempty(s1)
               Surfaceliferect=[640,385,600,550];
           else
               Surfaceliferect=[641,276,600,550];
           end 
    end
     if Diameter==2
       SurfacelifeArealimit=[90,1000];
        if ~isempty(s1)
           Surfaceliferect=[674,309,600,550];
       else
           Surfaceliferect=[641,276,600,550];
        end
     end
    if Diameter==3
       SurfacelifeArealimit=[100,1000];
        if ~isempty(s1)
           Surfaceliferect=[685,295,600,540];
       else
           Surfaceliferect=[685,295,600,520];
        end
 
    end
    if Diameter==4
       SurfacelifeArealimit=[200,2000];
       if ~isempty(s1)
           Surfaceliferect=[676,267,600,540];
       else
           Surfaceliferect=[669,251,600,530];
        end
    end
     if Diameter==5
       SurfacelifeArealimit=[300,5000];
       if ~isempty(s1)
           Surfaceliferect=[652,269,600,550];
       else
           Surfaceliferect=[659,280,600,550];
        end
    end
end

%%
vn=172;%气泡上身视频提取个数
flag=exist('vn','var');
interval=[1,2,3,6,12,15,30,60,120,250,400];%%
interval=sort(interval,'descend');
Name='bubble parameter analysis';
switch Name
    case {'videoextractionfun'}
for i=1:length(FileDirDiameter)
    filefoldname{i}=FileDirDiameter(i).name;
    exp='(?<=FPS|FPS Test)\s\d+';
    [s,e,strmatch]=regexp(filefoldname{i},exp,'start','end','match');
    FPS=str2num(string(strmatch));
    Directorylist(i,:)={[destinationini,filefoldname{i}]};
    Directory=[Directorylist{i,:},'\'];
    videoextractionfun(Directory,Arealimit,videolimit)%%侧视时候，气泡有效视频帧提取程序
end
    case{'FPSdeterminationmeancode.m'}%%FPS敏感性分析程序
for k=1:length(interval)
for i=1:length(FileDirDiameter)
%     try
    filefoldname{i}=FileDirDiameter(i).name;
    exp='(?<=FPS|FPS Test)\s\d+';
    [s,e,strmatch]=regexp(filefoldname{i},exp,'start','end','match');
    FPS=str2num(string(strmatch));
    FPS=FPS/interval(k);
    Directorylist(i,:)={[destinationini,filefoldname{i}]};
    Directory=[Directorylist{i,:},'\'];
    cd(Directory)
    mp4num=length(dir('*.mp4'));
    if mp4num==0
       mp4num=length(dir('*.avi'));
    end
   if flag==0
    vn=mp4num; 
   end
    FPSprefun(vn,Directory,Magnification,Arealimit,FPS,interval(k));
%     clearvarsfun(vn,Directory)
    FPSpostfun(vn,Directory,FPS);
    run(Name)
end
end
    case {'Bubblenumberdetermintationmeancode.m'}%%侧面拍摄气泡数量确定程序
for j=1:length(vn_candidate)
    vn=vn_candidate(j);
    flag=exist('vn','var');
    pack
for i=1:length(FileDirDiameter)
%     try
    filefoldname{i}=FileDirDiameter(i).name;
    
    exp='(?<=FPS|FPS Test)\s\d+';
    [s,e,strmatch]=regexp(filefoldname{i},exp,'start','end','match');
    FPS=str2num(string(strmatch));
    Directorylist(i,:)={[destinationini,filefoldname{i}]};
    Directory=[Directorylist{i,:},'\'];
    cd(Directory)
    mp4num=length(dir('*.mp4'));
    if mp4num==0
       mp4num=length(dir('*.avi'));
    end
   if flag==0
    vn=mp4num; 
   end
run(Name)
end
end
    case {'bubble parameter analysis'}%气泡侧面运动特性分析程序
for i=1:length(FileDirDiameter)
%     try
    filefoldname{i}=FileDirDiameter(i).name;
    
    exp='(?<=FPS|FPS Test)\s\d+';
    [s,e,strmatch]=regexp(filefoldname{i},exp,'start','end','match');
    FPS=str2num(string(strmatch));
    Directorylist(i,:)={[destinationini,filefoldname{i}]};
    Directory=[Directorylist{i,:},'\'];
    cd(Directory)
    mp4num=length(dir('*.mp4'));
    if mp4num==0
       mp4num=length(dir('*.avi'));
    end
 
    vn=mp4num; %现在终于可以用中文了
   
%    videoextractionfun(Directory,Arealimit,videolimit);
%  FPSprefun1(vn,Directory,Magnification,Arealimit,FPS);
   FPSpostfun(vn,Directory,FPS);
end
    case{'surfacelifeextraction'}%%提取含有目标的有效视频
for i=1:length(FileDirDiameter)
    %     try
    filefoldname{i}=FileDirDiameter(i).name;
    exp='(?<=FPS|FPS Test)\s\d+';
    [s,e,strmatch]=regexp(filefoldname{i},exp,'start','end','match');
    FPS=str2num(string(strmatch));
    Directorylist(i,:)={[destinationini,filefoldname{i}]};
    Directory=[Directorylist{i,:},'\'];
   if length(dir([Directory,'*.mp4']))>videolimit||length(dir([Directory,'*.avi']))>videolimit
    disp('video number is satisfied')
    return
   end
Directoryini=[Directory,'initialvideo\'];%特别注意这里的initialvideo指的是下面的文件
cd(Directoryini)
if exist('n_ini.txt','file')
    n_ini=textread('n_ini.txt');
else
    n_ini=1;
end
%% 对于所有的文件进行排序按照数字大小进行排序

filevideo=dir([Directoryini,'*.mp4']);%%need to change here
if length(filevideo)==0
    filevideo=dir([Directoryini,'*.avi']);
end
vftemp=length(filevideo);
n_end=vftemp;

[~,ind]=sort([filevideo(:).datenum],'ascend');
filevideo=filevideo(ind);

%% try to read the background
if exist([Directory,'background.jpg'])
    background=imread('background.jpg');
    background=imcrop(background,Surfaceliferect);
    imwrite(background,[Directory,'background.jpg']);
    [ylimit,xlimit]=size(background);
end
%% 识别跟我们要求的视频
for n=n_ini:n_end %n 代表有多少个拍摄的大视频
    %%
    fname=VideoReader([Directoryini,filevideo(n).name]);
    num=fname.NumberOfFrames;
    num_avi=length(dir([Directory,'*.avi']));
    if length(dir([Directory,'*.mp4']))>videolimit||length(dir([Directory,'*.avi']))>videolimit
        disp('video number is satisfied')
        return
    end
    clear videodetection Index   
%     h=waitbar(0,'Please wait the whole process');
    j=1;
    Index.ini=1;
    Index.end=num;
    for ni=1:num
        gray=read(fname,ni);
         gray=imcrop(gray,Surfaceliferect);
        [a,b,c]=size(gray);
        if c>1
            gray=rgb2gray(gray);
        end
%         imwrite(gray,[Directory,'gray.jpg']);
        if (ni==1)&&(~exist([Directory,'background.jpg']))
            background=read(fname,ni);
             if c>1
            background=rgb2gray(background);
             end
             background=imcrop(background,Surfaceliferect);
%             imwrite(background,[Directory,'background.jpg']);
            [ylimit,xlimit]=size(background);
        end
        if ni==1
            background=read(fname,ni);
            [a,b,c]=size(background);
        if c>1
            background=rgb2gray(background);
        end
           background=imcrop(background,Surfaceliferect);
        end
        
        %%
        graypicture=abs(gray-background);%
%          imwrite(graypicture,[Directory,'graypicture.jpg']);
        highlimit=150/255;% need to modify this number according to the light source intensity
        lowlimit=30/255;%nee to modify this number according to the light source intensity
        graycontrast=imadjust(graypicture,[lowlimit,highlimit],[0,1]);
%          imwrite(graycontrast,[Directory,'graycontrast.jpg']);
%%
flag=3;
switch flag
    case {1} 
        se90 = strel('line',3,90);
        se0 = strel('line',3,0); 
        graybinary=(graycontrast>75);
        graydilate=imdilate(graybinary,[se0 se90]);
        grayfill=imfill(graydilate,'holes');
        grayerode=imerode(grayfill,se0);
        grayerode=imerode(grayerode,se90);
    case {2}
        [~,threshold]=edge(graycontrast,'sobel');
        fudgefactor=0.5;
        graybinary=edge(graycontrast,'sobel',threshold*fudgefactor);
        se90 = strel('line',2,90);
        se0 = strel('line',2,0);
      graydilate=imdilate(graybinary,[se0 se90]);
        grayfill=imfill(graydilate,'holes');
          grayerode=imerode(grayfill,se0);
        grayerode=imerode(grayerode,se90);  
    case {3}
        graybinary=edge(graycontrast,'canny');
        se90 = strel('line',6,90);
        se0 = strel('line',6,0);
        graydilate=imdilate(graybinary,[se0 se90]);
        grayfill=imfill(graydilate,'holes');
        grayerode=imerode(grayfill,se0);
        grayerode=imerode(grayerode,se90);   
    case {4}
        graybinary=edge(graycontrast,'roberts');
        se90 = strel('line',2,90);
        se0 = strel('line',2,0);
        graydilate=imdilate(graybinary,[se0 se90]);
        grayfill=imfill(graydilate,'holes');
         grayerode=imerode(grayfill,se0);
         grayerode=imerode(grayerode,se90);
        
    case {5}
        graybinary=edge(graycontrast,'prewitt');
        se90 = strel('line',2,90);
        se0 = strel('line',2,0);
        graydilate=imdilate(graybinary,[se0 se90]);
        grayfill=imfill(graydilate,'holes');
        grayerode=imerode(grayfill,se0);
        grayerode=imerode(grayerode,se90);
 
    otherwise
        graybinary=imbinarize(graycontrast);
        grayerode=graybinary;
end

grayfill=grayerode;
% imshow(grayfill);
grayfilt=bwareafilt(grayfill,SurfacelifeArealimit);
grayfilt=imclearborder(grayfilt);
%% 处理Diameter03所需要的程序段
if(0)
        graybinary=graycontrast>75;%use imageproplus to get the best limit
% %          imwrite(graybinary,[Directory,'graybinary.jpg']);
        se=strel('disk',8);
        grayclose=imclose(graybinary,se);
        grayfill=imfill(grayclose,'holes');
%         imwrite(grayfill,[Directory,'grayfill.jpg']);
        grayfilt=bwareafilt(grayfill,SurfacelifeArealimit);
         grayfilt=imclearborder(grayfilt);
end
         %%

        frame=bwconncomp(grayfilt,8);
        % object=regionprops(frame,'all')
        videodetection(ni).number=frame.NumObjects;
        
    end
    j=1;
    Dvalue=6;% to produce Dvalue consecutive seeds which is used to detect how many images that have the same number of objects 
    for ki=1:length(videodetection)
        if ki<Dvalue+1||ki>(length(videodetection)-Dvalue-1)
            continue
        end
        if all([videodetection(ki-1).number]==0)&&all([videodetection(ki:ki+Dvalue).number]>0)
            Index(j).ini=ki-Dvalue;
        end
        if all([videodetection(ki-Dvalue:ki-1).number]>0)&&all([videodetection(ki:ki+Dvalue).number]==0)
            Index(j).end=ki+Dvalue;
            if Index(j).end-Index(j).ini<30
                Index(j)=[];
                j=j-1;
            end
            j=j+1;
        end
        
    end
    %% for some situation their is two different rising bubble in the Index_ini and Index_end which is not the ideal situation we want to have
   
  len=length(Index);%index is used to record the initial and final index correspond to bubble emergers and departure
    k=1;
    while(1)
        if len<=0
            break
        end
        if k>len
            break
        end
        Index_ini=Index(k).ini;
        Index_end=Index(k).end;
        if isempty(Index_end)
            Index(k)=[];
            len=len-1;
           if k>len
            break
           end
            Index_ini=Index(k).ini;
            Index_end=Index(k).end;
            
        end
        if k==1
            if isempty(Index_ini)
                Index(k).ini=1;
            end
    
        else
        if isempty(Index_ini) 
            Index(k)=[];
            len=len-1;
             if k>len
            break
             end
            Index_ini=Index(k).ini;
            Index_end=Index(k).end;
        end      
        end
        
        
        
        if ~isempty(Index_ini)&&~isempty(Index_end)%这里可能需要更因为有可能导致第一个视频过长，包含不想要的视频
         if (Index_ini==1)&&(Index_end==num)
        Index(k)=[];
        len=len-1;
         if k>len
            break
        end
        Index_ini=Index(k).ini;
        Index_end=Index(k).end;
         end
        end
         
        if ~isempty(Index(k).ini)&&~isempty(Index(k).end)
            k=k+1;
        end
    end
    %%
    k=1;
    del=[];
 for km=1:length(Index)
       if k>length(Index)
            break
       end
       if isempty(Index)
           break
       end
        Index_ini=Index(km).ini;
        Index_end=Index(km).end;  
        for kn=Index_ini+Dvalue:Index_end-Dvalue-1
            if kn+Dvalue>Index_end-Dvalue-1
                continue
            end
            if all([videodetection(kn:kn+Dvalue).number]==0)%要记录的视频中连续Dvalue张是0，就删除，中间不能有0值
                del(k)=km;
                k=k+1;
            end
        end
 end

 Index(del)=[];
    %%
    for m=1:length(Index)
        video=VideoWriter([Directory,num2str(m+num_avi),'.avi']);
        open(video);
        Index_i=Index(m).ini;
        Index_e=Index(m).end;
        if isempty(Index_i)||isempty(Index_e)
            break
        end
        for o=Index_i:Index_e
             saveframe=read(fname,o);
             saveframe=imcrop(saveframe,Surfaceliferect);
            writeVideo(video,saveframe);
        end
    end
    
    close(video)
    n_ini=n_ini+1;%pay more attention to this value, it represent which video it deal with.
    save([Directoryini,'n_ini.txt',],'n_ini','-ascii');
    
end
 
end
disp('finish')

case {'surfacelifeanalysis'}
%

 end





function FPSprefun1(vn,Directory,Magnification,Arealimit,FPS,interval)
cd(Directory)
filevideo=dir([Directory,'\*.mp4']);
if length(filevideo)==0
    filevideo=dir([Directory,'\*avi']);
end


%% 对于所有的文件进行排序，为什么这么难用
vftemp=length(filevideo);
for i=1:vftemp
   A(i)=str2num(filevideo(i).name(1:end-4)); 
end
[A,o]=sort(A);
for i=1:vftemp
    file(i)=filevideo(o(i));
end
filevideo=file;
if exist('background.jpg','file')
background=imread('background.jpg');
[ylimit,xlimit]=size(background);
else
    fname=VideoReader([Directory,filevideo(1).name]);
    background=read(fname,1);
    [a,b,c]=size(background);
    if c>1
    background=rgb2gray(background);
    end
    [ylimit,xlimit]=size(background);
    imwrite(background,[Directory,'background.jpg']);
end

for n=1:vn
   clear frameboundary fname gray dz dx object state detection vz vx delta_z delta_x vmagnititude 
%     try
fname=VideoReader([Directory,filevideo(n).name]);
num=fname.NumberOfFrames;
if nargin<6
    interval=1;
end
num=floor(num/interval);
%%
if  ~isempty(dir(['Bubble Result',fname.Name(1:end-4),'\','Bubble',' ','FPS',num2str(FPS),' ',num2str(n),'.mat']))
    continue  %detect if there is any destination folder, if exist, jump this cycle
end
%%
mkdir(['Bubble Result',fname.Name(1:end-4)])
h=waitbar(0,'Please wait the whole process');

h1=figure(1);
set(h1,'visible','off')
imshow(background);

% h2=figure(2);
% set(h2,'visible','on')
% imshow(background)
for i=1:num
%      if i==4
%         pause
%     end
    i_index=(i-1)*interval+1;
    gray=read(fname,i_index);
    if i_index==1
    backgroundshow=background;
    end
    [a,b,c]=size(gray);
    if c>1
    gray=rgb2gray(gray);
    end
    
    graypicture=abs(background-gray);%因为背景是亮的，所以是背景减去图片
%     imshow(graypicture)
    highlimit=50/255;
    lowlimit=10/255;
    graycontrast=imadjust(graypicture,[lowlimit,highlimit],[0,1]); 
%     imshow(graycontrast);
flag=1;
switch flag
    case {1} 
       
        se90 = strel('line',2,90);
        se0 = strel('line',2,0); 
        graybinary=(graycontrast>50);
        graydilate=imdilate(graybinary,[se0 se90]);
        grayfill=imfill(graydilate,'holes');
        grayerode=imerode(grayfill,se0);
        grayerode=imerode(grayerode,se90);
    case {2}
        [~,threshold]=edge(graycontrast,'sobel');
        fudgefactor=0.5;
        graybinary=edge(graycontrast,'sobel',threshold*fudgefactor);
        se90 = strel('line',2,90);
        se0 = strel('line',2,0);
      graydilate=imdilate(graybinary,[se0 se90]);
        grayfill=imfill(graydilate,'holes');
          grayerode=imerode(grayfill,se0);
        grayerode=imerode(grayerode,se90);  
    case {3}
        graybinary=edge(graycontrast,'canny');
        se90 = strel('line',2,90);
        se0 = strel('line',2,0);
        graydilate=imdilate(graybinary,[se0 se90]);
        grayfill=imfill(graydilate,'holes');
        grayerode=imerode(grayfill,se0);
        grayerode=imerode(grayerode,se90);   
    case {4}
        graybinary=edge(graycontrast,'roberts');
        se90 = strel('line',2,90);
        se0 = strel('line',2,0);
        graydilate=imdilate(graybinary,[se0 se90]);
        grayfill=imfill(graydilate,'holes');
         grayerode=imerode(grayfill,se0);
         grayerode=imerode(grayerode,se90);
        
    case {5}
        graybinary=edge(graycontrast,'prewitt');
        se90 = strel('line',2,90);
        se0 = strel('line',2,0);
        graydilate=imdilate(graybinary,[se0 se90]);
        grayfill=imfill(graydilate,'holes');
        grayerode=imerode(grayfill,se0);
        grayerode=imerode(grayerode,se90);
 
    otherwise
        graybinary=imbinarize(graycontrast);
        grayerode=graybinary;
end

     
     grayfilt=bwareafilt(grayerode,Arealimit);
%      h1=figure(2);
%      imshow(grayfilt)
%      delete(h1)

     frame=bwconncomp(grayfilt,8);
     if frame.NumObjects>1
     grayfilt=imclearborder(grayfilt);% clear the bubble that connect to the boundary which is very important to the postprocess
     frame=bwconncomp(grayfilt,8);
     end
%%

%% if we find the moving bubble touch the boudary, we will elimate these information

frameboundary=bwboundaries(grayfilt,'noholes');
if ~isempty(frameboundary)
bubbleboundary=frameboundary{1};%1st column is y coordinitor, 2nd column is x coordinitor
bubble_xmax=max(bubbleboundary(:,2));
bubble_xmin=min(bubbleboundary(:,2));
bubble_ymax=max(bubbleboundary(:,1));
bubble_ymin=min(bubbleboundary(:,1));
dz(i)=max(bubbleboundary(:,1))-min(bubbleboundary(:,1));
dx(i)=max(bubbleboundary(:,2))-min(bubbleboundary(:,2));
end
% end

state=regionprops(frame,'Area','Perimeter','MajorAxis','MinorAxis','Centroid','Orientation');
framepixelindex=frame.PixelIdxList;
[ms,~]=size(state);
%sometimes uwilling small bubbles will produce concommitant with the rising
%target bubble, therefore we need to think out a way to exclude this
%condition
%%
if (~ms||bubble_ymax==ylimit) 
object(i).Numobjects=0;
object(i).name=[num2str(i),'.jpg'];
continue
end
%%
if bubble_ymin==1% if the bubble touch the surface, all the information except the name of the frame will be clear 
    % note that when the bubble reaches the surface ,the boundary
    % coordinate will be 1 not ylimit
    for k=i:num
        object(k).name=[num2str(i),'.jpg'];
        object(k).Numobjects=0;
    end
    break%如果气泡已经解除到液面，就不需要再进行计算了
end
%%
% object(i).PixelIdxList=temp_object.PixelIdxList;
%continue if the Bubble Frameant image produced from A1 original image  subtracted by
%background image consists of wholy zero element, it could not be saved to
%current Directory in that the code will continue, or in another
%word,jumped to another cycle.
ns=size((fieldnames(state)),1);
s=zeros(ms,ns+3);
s(:,1:7)=[cat(1,state.Area),cat(1,state.Perimeter),cat(1,state.MajorAxisLength),cat(1,state.MinorAxisLength),cat(1, state.Centroid),cat(1,state.Orientation)];

s(:,ns+2)=4*pi*s(:,1)./(s(:,2).^2);%caculate the roundness of the bubble
s(:,ns+3)=s(:,3)./s(:,4);%此前的程序可能有错误
[sm,sn]=size(s);
if sm>1
    [maxs,maxi]=max(s(:,1));
    s_temp=s(maxi,:);
    s=s_temp;
end

object(i).name=[num2str(i),'.jpg'];
object(i).Numobjects=length(s(1,1));
object(i).Area=s(:,1);
object(i).Perimeter=s(:,2);
object(i).MajorAxisLength=s(:,3);
object(i).MinorAxisLength=s(:,4);
object(i).Centroidx=s(:,5);
object(i).Centroidy=s(:,6);
object(i).Orientation=s(:,7);
object(i).roudness=s(:,ns+2);
object(i).ratio=s(:,ns+3);
object(i).diameter=2*sqrt(object(i).Area/pi);
hold on 
   
    
if object(i).Numobjects~=0&&mod(i,8)==0 
    indexshow=(find(grayfilt(:)==1));
    backgroundshow(indexshow)=gray(indexshow);
%    imshow(backgroundshow);
   plot(object(i).Centroidx,object(i).Centroidy,'y*')
    hold on
end
   waitbar(i/num,h)
    if object(i).Numobjects==0
        continue
    end
    if n==54||n==53||n==52||n==96||n<6
        mkdir(['Bubble Frame',' ','FPS',num2str(FPS),' ',fname.Name(1:end-4)])
        imwrite(grayfilt,[Directory,'Bubble Frame',' ','FPS',num2str(FPS),' ',fname.Name(1:end-4),'\',object(i).name])
    end
    
end
close(h)

print('-dtiff',['Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPS),' ','trajectory.tif']);%
delete(h);
hold off
delete(h1)
h2=figure(2);
set(h2,'visible','off')
imshow(backgroundshow);
print('-dtiff',['Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPS),' ','trajectory1.tif']);%
delete(h2)
hold off
index=1;
%% find the frame when bubble occur or keep away from the screenshot 
detection(index).initial=1;
detection(index).end=num;
detection(index).delt=detection(index).end-detection(index).initial;
r_detection=length(detection);

%% 用来求得计算速度的首帧和末帧分别时多少
for j=1:num-1
    if object(1).Numobjects>0&&object(num).Numobjects>0
         break%%这里用来发现如果首帧和末帧气泡个数＞0就将计算的首末帧限制在1：num之间
    end
    if (object(j).Numobjects==0)&&(object(j+1).Numobjects>0)
    detection(index).initial=j+1;
    end
    if object(j).Numobjects>0&&object(j+1).Numobjects==0
    detection(index).end=j;
    detection(index).delt=detection(index).end-detection(index).initial;
    index=index+1;
    end
end
%% caculate the real velocity of bubble motion
dd=[detection.delt];
[sm,sn]=size(dd);
if sm>1||sn>1
%     pause
    [ni,nj]=max(dd);
    for ii=1:length(detection)
        if ii==nj
            continue
        end
        field=fieldnames(object);
        mf=max(size(field));
        for j=detection(ii).initial:detection(ii).end
        for i=1:mf
            if (strcmp(field{i},'name')||strcmp(field{i},'Numobjects'))
                continue
            end
            eval(['object(j).',field{i},'=[];'])
        end
        end
%     detection(ii)=[];
    end
    detection1=detection(nj);
    clear detection
    detection=detection1;
end

temp=detection.delt+1;
ini=detection.initial;
final=detection.end;
for t=ini:final-1
vz(t)= -(object(t+1).Centroidy-object(t).Centroidy)*Magnification/(1/FPS);%使用500FPS帧率拍摄需要确定到底放大倍数是多少,这里的0.1代表放大倍数mm/pixel
vx(t)= (object(t+1).Centroidx-object(t).Centroidx)*Magnification/(1/FPS);%计算速度是多少
delta_z(t)=(object(t+1).Centroidy-object(t).Centroidy)*Magnification;%计算当前气泡的距离是多少
delta_x(t)=(object(t+1).Centroidx-object(t).Centroidx)*Magnification;%0.1代表0.1mm/pixel
end
%% Make sure the robustness of this code
try

if final<num-1
vz(:,final:end)=[];
vx(:,final:end)=[];
delta_x(:,final:end)=[];
end
if ini>1
vz(:,1:ini-1)=[];
vx(:,1:ini-1)=[];
delta_z(:,1:ini-1)=[];
end
catch
end
%% get the statistic of the bubbles
[a,b,c]=size(background);
cz=(a-cat(1,object.Centroidy))*Magnification;
cx=(cat(1,object.Centroidx)-b/2)*Magnification;
majoraxis=cat(1,object.MajorAxisLength)*Magnification;
minoraxis=cat(1,object.MinorAxisLength)*Magnification;
area=cat(1,object.Area)*Magnification^2;
diameter_area=sqrt(4.*area./pi);
diameter_axis=(majoraxis.^2.*minoraxis).^(1/3);
aspectratio=minoraxis./majoraxis;

x=linspace(1/FPS,(final-ini)/FPS,final-ini);
vmagnititude=sqrt(vz.^2+vx.^2);
cz(1,:)=[];cx(1,:)=[];aspectratio(1,:)=[];diameter_axis(1,:)=[];diameter_area(1,:)=[];majoraxis(1,:)=[];minoraxis(1,:)=[];area(1,:)=[];
velocity_ins=[cz,cx,vz',aspectratio,diameter_axis,diameter_area,vx',vmagnititude',majoraxis,minoraxis,area];
[mi,ni]=size(velocity_ins);
velocity_ins_cell=mat2cell(velocity_ins,ones(mi,1),ones(ni,1));
velocity_ins_title={'cz','cx','vz_ins','aspectratio','diameter_volume','diameter_area','vx_ins','vm_ins','majoraxis','minoraxis','area'};
velocity_ins_result=[velocity_ins_title;velocity_ins_cell];
xlswrite([Directory,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPS),' ','velocity_ins.xls'],velocity_ins_result);


vx_ave=mean(vx);  
vm_ave=mean(vmagnititude);
vz_ave=mean(vz); 
majoraxis_ave=mean(majoraxis);
minoraxis_ave=mean(minoraxis);
area_ave=mean(area);
diameter_area_ave=mean(diameter_area);
diameter_axis_ave=mean(diameter_axis);
aspectratio_ave=mean(aspectratio);

std_majoraxis=std(majoraxis);
std_minoraxis=std(minoraxis);
std_area=std(area);
std_diameter_area=std(diameter_area);
std_diameter_axis=std(diameter_axis);
std_aspectratio=std(aspectratio);
std_vx=std(vx);
std_vz=std(vz);
std_vm=std(vmagnititude);
velocity_mean=[vz_ave,aspectratio_ave,diameter_axis_ave,diameter_area_ave,vz_ave,vm_ave,majoraxis_ave,minoraxis_ave,area_ave];
std_v=[std_vz,std_aspectratio,std_diameter_axis,std_diameter_area,std_vz,std_vm,std_majoraxis,std_minoraxis,std_area];
velocity_mean=[velocity_mean,std_v];
[mv,nv]=size(velocity_mean);
velocity_mean_cell=mat2cell(velocity_mean,ones(mv,1),ones(nv,1));
velocity_mean_title={'vz_ave','aspectratio_ave','diameter_axis_ave','diameter_area_ave','vz_ave','vm_ave','majoraxis_ave','minoraxis_ave','area_ave'...
    'std_vz','std_aspectratio','std_diameter_axis','std_diameter_area','std_vz','std_vm','std_majoraxis','std_minoraxis','std_area'};
velocity_mean_result=[velocity_mean_title;velocity_mean_cell];

save([Directory,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPS),' ','velocity_mean.txt'],'velocity_mean','-ascii')
xlswrite([Directory,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPS),' ','velocity_mean.xls'],velocity_mean_result);
save([Directory,'Bubble Result',fname.Name(1:end-4),'\','Bubble',' ','FPS',num2str(FPS),' ',fname.Name(1:end-4),'.mat'])

end
end

%%
function FPSpostfun(vn,Directory,FPS)
FPStemp=FPS;
vntemp=vn;
Directorytemp=Directory;
cd(Directorytemp)
filelist=dir([Directorytemp,'Bubble Result*']);%读取所的


%% 待会儿需要删除这一段
vn=length(filelist);
vntemp=length(filelist);    
%% set the order of the files
for i=1:vntemp
   filelist_name=filelist(i).name;
   sn(i)=str2num(filelist_name(isstrprop(filelist_name,'digit'))); % if the number of files is over 100,change end-1 to end-2
end
[sn,n]=sort(sn);
for i=1:vntemp
    tempfile(i)=filelist(n(i));
end
filelist=tempfile;
%% plot the mean velocity vs time
for ti=1:vntemp
  
    cd([Directorytemp, filelist(ti).name])
    fmat=dir(['Bubble',' ','FPS',num2str(FPStemp),' ',num2str(ti),'.mat']);
    fmat=dir('*.mat');

    load(fmat.name);  
    %% velocity information

    x=linspace(1/FPS,(final-ini)/FPS,final-ini);
    vmagnititude=sqrt(vz.^2+vx.^2);
    vx_ave=mean(vx);%vx means the velocity of different frames
    vm_ave=mean(vmagnititude);
    vz_ave=mean(vz);%vz means the velocity of different frames
    
    Velocity=[vx_ave;vz_ave;vm_ave];
    std_vx=std(vx);
    std_vz=std(vz);
    std_vm=std(vmagnititude);
    std_v=[std_vx,std_vz,std_vm];
%% axis information
    majoraxis=cat(1,object.MajorAxisLength)*Magnification;
    minoraxis=cat(1,object.MinorAxisLength)*Magnification;
    area=cat(1,object.Area)*Magnification^2;
    diameter_area=sqrt(4.*area./pi);
    diameter_axis=(majoraxis.^2.*minoraxis).^(1/3);
    aspectratio=minoraxis./majoraxis;
    %%
    majoraxis_ave=mean(majoraxis);
    minoraxis_ave=mean(minoraxis);
    area_ave=mean(area);
    diameter_area_ave=mean(diameter_area);
    diameter_axis_ave=mean(diameter_axis);
    aspectratio_ave=mean(aspectratio);
    %%
    std_majoraxis=std(majoraxis);
    std_minoraxis=std(minoraxis);
    std_area=std(area);
    std_diameter_area=std(diameter_area);
    std_diameter_axis=std(diameter_axis);
    std_aspectratio=std(aspectratio);
    std_parameter=[std_majoraxis std_minoraxis std_area std_diameter_area std_diameter_axis std_aspectratio];% the order to writer into a txt file
    %is std_majoraxis std_minoraxis std_area std_diameter_area
    %std_diameter_axis
    ave_parameter=[majoraxis_ave minoraxis_ave area_ave,diameter_area_ave diameter_axis_ave aspectratio_ave];
   %% statistic information of velocity, major axis, minor axis, area, equiavelent diameter. summarize it into velocity_bubbles structure
    velocity_bubbles(ti).vxmean=vx_ave;%%A(1,1);
    velocity_bubbles(ti).vzmean=vz_ave;%%A(1,2);
    velocity_bubbles(ti).vmmean=vm_ave;%%A(1,3);
    diameter_bubbles(ti).areamean=diameter_area_ave;
    diameter_bubbles(ti).axismean=diameter_axis_ave;
    axis_bubbles(ti).majoraxismean=majoraxis_ave;
    axis_bubbles(ti).minoraxismean=minoraxis_ave;
    area_bubbles(ti).areamean=area_ave;
    aspectratio_bubbles(ti).aspectratiomean=aspectratio_ave;
    deltat(ti,:)=(final-ini+1)./FPS;
      seed=1;
      clear h
   %%
    save([Directorytemp,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPStemp),' ','std parameter.txt'],'std_parameter','-ascii')
    save([Directorytemp,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPStemp),' ','ave parameter.txt'],'ave_parameter','-ascii')
    %% find the amplitude,cycle time,oscillation amplitude
    cxnew=cx-cx(1);
    cxmean=movmean(cx,1);
    cxmean=smooth(cxnew);
    cxmean=cxmean-cxmean(1);
    t=((1:1:final-ini)./FPS)';
    
    try
    [maxpks,maxlcs]=findpeaks(cxmean,t,'MinPeakDistance',0.1,'MinPeakHeight',2);
    [minpks,minlcs]=findpeaks(-cxmean,t,'MinPeakDistance',0.1,'MinPeakHeight',2);
    
    catch    disp([filelist(ti).name,' check']);
        continue
    end
    minpks=-minpks;
    [smaxpks,~]=size(maxpks);
    [sminpks,~]=size(minpks);
    if smaxpks>1
        cycleT=mean(diff(maxlcs));
    else
        if sminpks>1
            cycleT=mean(diff(minlcs));
        end
    end
    if (smaxpks==1&&sminpks==0)||(smaxpks==0&&sminpks==1)
%         disp([num2str(ti),' has a problem']);
    end
     if (smaxpks==2&&sminpks==1)
        cycleT=mean(diff(maxlcs));
    end
    if (smaxpks==1&&sminpks==2)
        cycleT=mean(diff(minlcs));
    end
    if smaxpks==sminpks
        cycleT=mean(abs(maxlcs-minlcs)*2);
    end    
    osctime=min(smaxpks,sminpks);
    if smaxpks>sminpks
   oscamp=maxpks(1:sminpks)-minpks(1:sminpks);     
    else 
   oscamp=maxpks(1:smaxpks)-minpks(1:smaxpks);
    end
    if size(oscamp,1)==0
       oscamp=0; 
    end
     if size(cycleT,1)==0
       cycleT=0; 
     end
  
    toscamp=sum(oscamp);
    %% find the wave length
      
   try
    [maxpksz,maxlcsz]=findpeaks(cxmean,cz,'MinPeakDistance',15,'MinPeakHeight',2);
    [minpksz,minlcsz]=findpeaks(-cxmean,cz,'MinPeakDistance',15,'MinPeakHeight',2);
    
    
   catch    disp([filelist(ti).name,' check']);
       continue
       
   end
    minpksz=-minpksz;
    [smaxpksz,~]=size(maxpksz);
    [sminpksz,~]=size(minpksz);
    if smaxpksz>1
        cycleL=mean(diff(maxlcsz));
    else
        if sminpksz>1
            cycleL=mean(diff(minlcsz));
        end
    end
    if (smaxpksz>0&&sminpksz==0)||(smaxpksz==0&&sminpksz>0)
%         disp([num2str(ti),' has a problem']);
    end
     if (smaxpksz==2&&sminpksz==1)
        cycleL=mean(diff(maxlcsz));
    end
    if (smaxpksz==1&&sminpksz==2)
        cycleL=mean(diff(minlcsz));
    end
    if smaxpksz==sminpksz
        cycleL=mean(abs(maxlcsz-minlcsz)*2);
    end 
    if size(cycleL,1)==0
       cycleL=0; 
    end
    oscillation(ti).cycleT=cycleT;
    oscillation(ti).cycleL=cycleL;
    oscillation(ti).toscamp=toscamp;
    oscillation(ti).osctime=osctime;
    oscillation(ti).oscamp=toscamp./osctime;
    %%
   if ~isempty(dir([Directorytemp,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPStemp),' ','vz.tif']))
    continue
    end


 %% plot cx vs t
     clf('reset');close all
    h=figure(seed);
     set(h,'visible','on')
    plot(x,cx-cx(1),'g-s','MarkerFaceColor','g','MarkerEdgeColor','g')
    % plot each frame velocity
    hold on
    plot(minlcs,minpks,'o')
    hold on 
    plot(maxlcs,maxpks,'o')
    title('instaneous X vs t')
    xlabel('time/s')
    ylabel('X/mm')
    legend('X','Minimum','Maxmum','Location','best')
    hold off
    print('-dtiff',[Directorytemp,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPStemp),' ','X vs t.tif']);
    saveas(gcf,['FPS',num2str(FPStemp),' ','X vs t.fig'])
    delete(h)
    seed=seed+1;
    %% plot cx vs cz
     clf('reset');close all
    h=figure(seed);
     set(h,'visible','on')
    plot(cz,cx-cx(1),'g-s','MarkerFaceColor','g','MarkerEdgeColor','g')
    % plot each frame velocity
    hold on
    plot(minlcsz,minpksz,'o')
    hold on 
    plot(maxlcsz,maxpksz,'o')
    title('instaneous X vs t')
    xlabel('Z/mm')
    ylabel('X/mm')
    legend('X','Minimum','Maxmum','Location','best')
    hold off
    print('-dtiff',[Directorytemp,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPStemp),' ','X vs Z.tif']);
    saveas(gcf,['FPS',num2str(FPStemp),' ','X vs Z.fig'])
    delete(h)
    seed=seed+1;
     if ~isempty(dir([Directorytemp,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPStemp),' ','vz.tif']))
    continue
    end
    
   %% start plotting minor axis and major axis together
    xs=1/FPS:1/FPS:(final-ini+1)/FPS;
    clf('reset');close all
    h=figure(seed);
     set(h,'visible','on')
    plot(xs,minoraxis,'k-s','MarkerFaceColor','k','MarkerEdgeColor','k')           % plot minor axis
    hold on
    plot(xs,majoraxis,'r-o','MarkerFaceColor','r','MarkerEdgeColor','r')
    hold on
    plot([1/FPS,(final-ini+1)/FPS],[minoraxis_ave,minoraxis_ave],'--k')   % Plot average minor axis
    hold on
    plot([1/FPS,(final-ini+1)/FPS],[majoraxis_ave,majoraxis_ave],'--r')
    title('major axis and minor axis vs time')
    xlabel('time/s')
    ylabel('major axis/mm')
    legend('minor axis','major axis','mean minor axis','mean major axis','Location','best')
    hold off
    print('-dtiff',[Directorytemp,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPStemp),' ','minor and major axis.tif']);
    saveas(gcf,['FPS',num2str(FPStemp),' ','minor and major axis.fig'])
    delete(h)
    seed=seed+1;
  
     %% comparison between major axis and minor axis
      clf('reset');close all
    h=figure(seed);
     set(h,'visible','on')
    scatter(minoraxis,majoraxis,'k','filled')           % plot minor axis
    xlabel('minor axis/mm')
    ylabel('major axis/mm')
    title('comparison of major axis and minor axis')
    hold off
    print('-dtiff',[Directorytemp,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPStemp),' ','major and minor axis 2.tif']);
     saveas(gcf,['FPS',num2str(FPStemp),' ','major and minor axis 2.fig'])
     delete(h)
   %% plot the area
    xs=1/FPS:1/FPS:(final-ini+1)/FPS;
    clf('reset');close all
    h=figure(seed);
     set(h,'visible','on')
    plot(xs,area,'g-s','MarkerFaceColor','g','MarkerEdgeColor','g')   % plot each frame velocity
    hold on
    plot([1/FPS,(final-ini+1)/FPS],[area_ave,area_ave],'--r')   % Plot average rising velocity
    title('bubble area')
    xlabel('time/s')
    ylabel('bubble area/mm^2')
    hold off
    print('-dtiff',[Directorytemp,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPStemp),' ','bubble area.tif']);
    saveas(gcf,['FPS',num2str(FPStemp),' ','bubble area.fig'])
    delete(h)
    seed=seed+1;
    %% plot double axis figure画双坐标图
     xs=1/FPS:1/FPS:(final-ini)/FPS;
    clf('reset');close all
    h=figure(seed);
     set(h,'visible','on')
    yyaxis left
    plot(xs,vz,'k-s','MarkerFaceColor','k','MarkerEdgeColor','k')           % plot vz
    hold on 
    plot(xs,vx,'r-o','MarkerFaceColor','r','MarkerEdgeColor','r')
    xlabel('time/s')
    ylabel('velocity/mm·s^{-1}')
    ax=gca;
    ax.YColor='k';
    yyaxis right
    aspectratio1=aspectratio;
    aspectratio1(end,:)=[];
    plot(xs,aspectratio1,'g-^')
    ylabel('aspect ratio');
    ax=gca;
    ax.YColor='g';
    legend('vz','vx','aspect ratio','Location','best')
    print('-dtiff',[Directorytemp,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPStemp),' ','velocity and aspect ration.tif']);
    saveas(gcf,['FPS',num2str(FPStemp),' ','velocity and aspect ratio.fig'])
    delete(h)
    seed=seed+1;
    
    
     xs=1/FPS:1/FPS:(final-ini)/FPS;
    clf('reset');close all
    h=figure(seed);
     set(h,'visible','on')
    yyaxis left
    plot(xs,vz,'k-s','MarkerFaceColor','k','MarkerEdgeColor','k')           % plot vz
    xlabel('time/s')
    ylabel('velocity/mm·s^{-1}')
    ax=gca;
    ax.YColor='k';
    yyaxis right
   aspectratio1=aspectratio;
    aspectratio1(end,:)=[];
    plot(xs,aspectratio1,'g-^')
    ylabel('aspect ratio');
    ax=gca;
    ax.YColor='g';
    legend('vz','aspect ratio','Location','best')
    print('-dtiff',[Directorytemp,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPStemp),' ','velocity and aspect ratio1.tif']);
    saveas(gcf,['FPS',num2str(FPStemp),' ','velocity and aspect ratio1.fig'])
    delete(h)
    seed=seed+1;
    
    
     xs=1/FPS:1/FPS:(final-ini)/FPS;
    clf('reset');close all
    h=figure(seed);
     set(h,'visible','on')
    yyaxis left
    plot(xs,vmagnititude,'k-s','MarkerFaceColor','k','MarkerEdgeColor','k')           % plot vz
    xlabel('time/s')
    ylabel('velocity/mm·s^{-1}')
    ax=gca;
    ax.YColor='k';
    yyaxis right
   aspectratio1=aspectratio;
    aspectratio1(end,:)=[];
    plot(xs,aspectratio1,'g-^')
    ylabel('aspect ratio');
    ax=gca;
    ax.YColor='g';
    legend('vm','aspect ratio','Location','best')
    print('-dtiff',[Directorytemp,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPStemp),' ','velocity and aspect ratio2.tif']);
    saveas(gcf,['FPS',num2str(FPStemp),' ','velocity and aspect ratio2.fig'])
    delete(h)
    seed=seed+1;
    
    
    
    xs=1/FPS:1/FPS:(final-ini)/FPS;
    clf('reset');close all
    h=figure(seed);
     set(h,'visible','on')
    scatter(diameter_axis,aspectratio,'k','filled')
    ylabel('aspectratio')
    xlabel('volume equivalent diameter/mm')
    print('-dtiff',[Directorytemp,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPStemp),' ','diameter axis and aspect ratio.tif']);
    saveas(gcf,['FPS',num2str(FPStemp),' ','diameter axis and aspect ratio.fig'])
    delete(h)
    seed=seed+1;
    
    
 
    %% bubble equivalent diameter comparison
     clf('reset');close all
    xs=1/FPS:1/FPS:(final-ini+1)/FPS;
    h=figure(seed);
     set(h,'visible','on')
    plot(xs,diameter_axis,'k-o','MarkerFaceColor','k','MarkerEdgeColor','k')   % plot each frame velocity
    hold on
    plot(xs,diameter_area,'r-s','MarkerFaceColor','r','MarkerEdgeColor','r')   % plot each frame velocity
    hold on
    plot([1/FPS,(final-ini+1)/FPS],[diameter_axis_ave, diameter_axis_ave],'k--')   % Plot average rising velocity
    hold on
    plot([1/FPS,(final-ini+1)/FPS],[diameter_area_ave, diameter_area_ave],'r--')   % Plot average rising velocity
    xlabel('time/s')
    ylabel('bubble equivalent diameter/mm')
    legend('equivalent volume diameter','equivalent area diameter','Location','best')
    hold off
    print('-dtiff',[Directorytemp,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPStemp),' ','bubble equivalent diameter comparison 2.tif']);
     saveas(gcf,['FPS',num2str(FPStemp),' ','bubble equivalent diameter comparison 2.fig'])
     delete(h)
    seed=seed+1;
    %% histogram of the area
    clf('reset');close all
    h=figure(seed);
     set(h,'visible','on')
    histogram(area,25)
    xlabel('area/mm^2')
    ylabel('frequency')
    title('area of bubble')
    print('-dtiff',[Directorytemp,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPStemp),' ','area of bubbles.tif'])
    saveas(gcf,['FPS',num2str(FPStemp),' ','area of bubbles.fig'])
    delete(h)
    seed=seed+1;   
    %% plot the diameter distribution of bubbles
    seed=seed+1;
    clf('reset');close all
    h=figure(seed);
     set(h,'visible','on')
    histogram(diameter_axis,25)
    hold on 
    histogram(diameter_area,25)
    legend('equivalent volume diameter', 'equivalent area diameter','Location','best')
    xlabel('equivalent diameter distribution/mm')
    ylabel('frequency')
    title('equivalent diameter')
    print('-dtiff',[Directorytemp,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPStemp),' ','equivalent diameter distribution comparison.tif'])
    saveas(gcf,['FPS',num2str(FPStemp),' ','equivalent diameter distribution comparison.fig'])
    delete(h)
    seed=seed+1;
    
    h=figure(seed);
set(h,'visible','on')
scatter(diameter_area,diameter_axis,'k','filled')
hold on
diametermax=max([max(diameter_area),max(diameter_axis)]);
diametermin=min([min(diameter_area),min(diameter_axis)]);
plot([diametermin,diametermax],[diametermin,diametermax],'--r')
xlabel('equivalent area diamter/mm')
ylabel('equivalent volume diameter/mm')
hold off
print('-dtiff',[Directorytemp,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPStemp),' ','equivalent diameter distribution comparison1.tif'])
    saveas(gcf,['FPS',num2str(FPStemp),' ','equivalent diameter distribution comparison1.fig'])
delete(h)
seed=seed+1;


    %% plot total velocity
    clf('reset');close all
    h=figure(seed);
     set(h,'visible','on')
    plot(x,vx,'k-s','MarkerFaceColor','k','MarkerEdgeColor','k')
    hold on
    plot(x,vz,'r-o','MarkerFaceColor','r','MarkerEdgeColor','r')
    hold on
    plot(x,vmagnititude,'g-d','MarkerFaceColor','g','MarkerEdgeColor','g')
    % plot each frame velocity
    hold on
    plot([0,(temp-2)/FPS],[vx_ave,vx_ave],'--k')   % Plot average lateral velocity
    hold on
    plot([0,(temp-2)/FPS],[vz_ave,vz_ave],'--r')
    hold on
    plot([0,(temp-2)/FPS],[vm_ave,vm_ave],'--g')
    hold on
    title('instaneous velocity of one bubble during rising process')
    xlabel('time/s')
    ylabel('velocity/mm·s^{-1}')
    legend('vx','vz','v magnititude','vx mean','vz mean','vm mean','Location','best')
    hold off
    print('-dtiff',[Directorytemp,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPStemp),' ','vtotal.tif']);
    saveas(gcf,['FPS',num2str(FPStemp),' ','vtotal.fig'])
    delete(h)
    seed=seed+1;
   
     %% comparison of vz and vx for one bubble
    clf('reset');close all
    h=figure(seed);
     set(h,'visible','on')
    scatter(vx,vz,'k','filled')
    xlabel('instanteneous v_x/mm·s^{-1}')
    ylabel('instanteneous v_z/mm·s^{-1}')
    title('comparison of instanteneous v_z and instanteneous v_x')
    
    hold off
    print('-dtiff',[Directorytemp,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPStemp),'comparison of instanteneous velocity_z and instanteneous velocity_x.tif']);
    saveas(gcf,['FPS',num2str(FPStemp),' ','comparison of instanteneous velocity_z and instanteneous velocity_x.fig'])
    delete(h)
    seed=seed+1;
    %%  histogram of the velocity
    clf('reset');close all
    h=figure(seed);
     set(h,'visible','on')
    histogram(vmagnititude,25)
    xlabel('v magnititude/mm·s^{-1}')
    ylabel('frequency')
    title('v magnititude of bubble')
    print('-dtiff',[Directorytemp,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPStemp),' ','v magnititude histogram.tif'])
    saveas(gcf,['FPS',num2str(FPStemp),' ','v magnititude histogram.fig'])
    delete(h)
    seed=seed+1;
    
    clf('reset');close all
    h=figure(seed);
     set(h,'visible','on')
    histogram(vz,25)
    xlabel('v_z/mm·s^{-1}')
    ylabel('frequency')
    title('v_z of bubble')
    print('-dtiff',[Directorytemp,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPStemp),' ','vz.tif'])
    saveas(gcf,['FPS',num2str(FPStemp),' ','vz.fig'])
    delete(h)
    seed=seed+1;
    
    clf('reset');close all
    h=figure(seed);
     set(h,'visible','on')
    histogram(vx,25)
    xlabel('v_x/mm·s^{-1}')
    ylabel('frequency')
    title('v_x of bubble')
    print('-dtiff',[Directorytemp,'Bubble Result',fname.Name(1:end-4),'\','FPS',num2str(FPStemp),' ','vx.tif'])
    saveas(gcf,['FPS',num2str(FPStemp),' ','vx.fig'])
    delete(h)
    seed=seed+1;
    
end
cd(Directorytemp)
  mkdir(['PostProcess',num2str(vntemp)])
%% plot the mean parameters over time
seed=1;
clear h
vxmean=cat(1,velocity_bubbles.vxmean);%time mean velocity of different bubble during rising processes
vzmean=cat(1,velocity_bubbles.vzmean);
vmmean=cat(1,velocity_bubbles.vmmean);
diameter_areamean=cat(1,diameter_bubbles.areamean);
diameter_axismean=cat(1,diameter_bubbles.axismean);
majoraxismean=cat(1,axis_bubbles.majoraxismean);
minoraxismean=cat(1,axis_bubbles.minoraxismean);
aspectratiomean=cat(1,aspectratio_bubbles.aspectratiomean);
cycleTmean=cat(1,oscillation.cycleT);
cycleLmean=cat(1,oscillation.cycleL);
toscampmean=cat(1,oscillation.toscamp);
oscampmean=cat(1,oscillation.oscamp);
oscampmean(isnan(oscampmean))=0;
osctimemean=cat(1,oscillation.osctime);
vmean_inf=[deltat,vzmean,aspectratiomean,diameter_axismean,diameter_areamean,vxmean,vmmean,majoraxismean,minoraxismean,cycleTmean,cycleLmean,toscampmean,oscampmean,osctimemean];
vmean_inf_title={'deltat','vzmean','aspectratiomean','diameter_volumemean','diameter_areamean','vxmean','vmmean','majoraxismean','minoraxismean','cycleTmean','cycleLmean','toscampmean','oscampmean','osctimemean'};
[mt,nt]=size(vmean_inf);
vmean_inf_cell=mat2cell(vmean_inf,ones(mt,1),ones(nt,1));
vmean_inf_result=[vmean_inf_title;vmean_inf_cell];
xlswrite([Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','vmean_inf.xls'],vmean_inf_result);

vxmean_bubbles=mean(vxmean);
vzmean_bubbles=mean(vzmean);
vmmean_bubbles=mean(vmmean);
std_vxmean_bubbles=std(vxmean);
std_vzmean_bubbles=std(vzmean);
std_vmmean_bubbles=std(vmmean);

diameter_areamean_bubbles=mean(diameter_areamean);
diameter_axismean_bubbles=mean(diameter_axismean);
majoraxismean_bubbles=mean(majoraxismean);
minoraxismean_bubbles=mean(minoraxismean);
aspectratiomean_bubbles=mean(aspectratiomean);
std_majoraxismean_bubbles=std(majoraxismean);
std_minoraxismean_bubbles=std(minoraxismean);
std_diameter_areamean_bubbles=std(diameter_areamean);
std_diameter_axismean_bubbles=std(diameter_axismean);
std_aspectratiomean_bubbles=std(aspectratiomean);
%% z
cycleTmean_bubbles=mean(cycleTmean(~cycleTmean==0));
cycleLmean_bubbles=mean(cycleLmean(~cycleLmean==0));
toscampmean_bubbles=mean(toscampmean(~toscampmean==0));
oscampmean_bubbles=mean(oscampmean(~oscampmean==0));
osctimemean_bubbles=mean(osctimemean(~osctimemean==0));

std_cycleTmean_bubbles=std(cycleTmean(~cycleTmean==0));
std_cycleLmean_bubbles=std(cycleLmean(~cycleLmean==0));
std_toscampmean_bubbles=std(toscampmean(~toscampmean==0));
std_oscampmean_bubbles=std(oscampmean(~oscampmean==0));
std_osctimemean_bubbles=std(osctimemean(~osctimemean==0));




xt=linspace(1,length(vxmean),length(vxmean));
vmean_infbubbles=[vzmean_bubbles,aspectratiomean_bubbles,diameter_axismean_bubbles,diameter_areamean_bubbles,vxmean_bubbles,vmmean_bubbles,majoraxismean_bubbles,minoraxismean_bubbles,...
                 cycleTmeanbubbles,cycleLmean_bubbles,toscampmean_bubbles,oscampmean_bubbles,osctimemean_bubbles,...
                std_vzmean_bubbles,std_aspectratiomean_bubbles,std_diameter_axismean_bubbles,std_diameter_areamean_bubbles,std_vxmean_bubbles,std_vmmean_bubbles,std_majoraxismean_bubbles,std_minoraxismean_bubbles,...
                std_cycleTmean_bubbles,std_cycleLmean_bubbles,std_toscampmean_bubbles,std_oscampmean_bubbles,std_osctimemean_bubbles];
[mt,nt]=size(vmean_infbubbles);
vmean_infbubbles_cell=mat2cell(vmean_infbubbles,ones(mt,1),ones(nt,1));
vmean_infbubbles_title={'vzmean_bubbles','aspectratiomean_bubbles','diameter_volumemean_bubbles','diameter_areamean_bubbles','vxmean_bubbles','vmmean_bubbles','majoraxismean_bubbles','minoraxismean_bubbles',...
                 'cycleTmean_bubbles','cycleTmean_bubbles','cycleLmean_bubbles','toscampmean_bubbles','oscampmean_bubbles','osctimemean_bubbles',...
               'std_vzmean_bubbles','std_aspectratiomean_bubbles','std_diameter_volumemean_bubbles','std_diameter_areamean_bubbles','std_vxmean_bubbles','std_vmmean_bubbles','std_majoraxismean_bubbles','std_minoraxismean_bubbles',...
               'std_cycleTmean_bubbles','std_cycleLmean_bubbles','std_toscampmean_bubbles','std_oscampmean_bubbles','std_osctimemean_bubbles'};
vmean_infbubbles_result=[vmean_infbubbles_title;vmean_infbubbles_cell];
xlswrite([Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','vmean_infbubbles.xls'],vmean_infbubbles_result);
Std_vemean_bubbles=[std_vxmean_bubbles,std_vzmean_bubbles,std_vmmean_bubbles];
Std_axismean_bubbles=[std_majoraxismean_bubbles,std_minoraxismean_bubbles,std_diameter_areamean_bubbles,std_diameter_axismean_bubbles,std_aspectratiomean_bubbles];
Vmean_bubbles=[vxmean_bubbles, vzmean_bubbles,vmmean_bubbles];
Axismean_bubbles=[majoraxismean_bubbles,minoraxismean_bubbles,diameter_areamean_bubbles,diameter_axismean_bubbles,aspectratiomean_bubbles];

save([Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','Vemean_bubbles.txt'],'Vmean_bubbles','-ascii')
save([Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','Axismean_bubbles.txt'],'Axismean_bubbles','-ascii')
save([Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','Std_vemean_bubbles.txt'],'Std_vemean_bubbles','-ascii')
save([Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','Std_axismean_bubbles.txt'],'Std_axismean_bubbles','-ascii')
 if isempty(dir([Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),'velocity_z distribution of different bubbles'])) 
     
%% plot average bubble rising velocity over different rising processes
clf('reset');close all
h=figure(seed);
 set(h,'visible','on')
plot(xt,vxmean,'k-s','MarkerFaceColor','k','MarkerEdgeColor','k')    
 hold on 
plot(xt,vzmean,'r-o','MarkerFaceColor','r','MarkerEdgeColor','r')
    hold on
plot(xt,vmmean,'g-d','MarkerFaceColor','g','MarkerEdgeColor','g')    
hold on 
plot([1,vntemp],[vxmean_bubbles,vxmean_bubbles],'--k')
hold on% Plot average lateral velocity
plot([1,vntemp],[vzmean_bubbles,vzmean_bubbles],'--r') 
hold on
plot([1,vntemp],[vmmean_bubbles,vmmean_bubbles],'--g') 
xlabel('n_{th} bubble')
ylabel('velocity/mm·s^{-1}')
title('mean velocity of different bubbles')
legend('vx mean','vz mean','vm mean','Location','best')
print('-dtiff',[Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','velocity of different bubbles.tif'])
saveas(gcf,[Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','velocity of different bubbles.fig'])
hold off 
delete(h)
seed=seed+1;
%% plot average bubble majoraxismean and monoraxis mean
clf('reset');close all
h=figure(seed);
 set(h,'visible','on')
plot(xt,majoraxismean,'k-s','MarkerFaceColor','k','MarkerEdgeColor','k')    
 hold on 
plot(xt,minoraxismean,'r-o','MarkerFaceColor','r','MarkerEdgeColor','r')  
hold on 
plot([1,vntemp],[majoraxismean_bubbles,majoraxismean_bubbles],'--k')
hold on% Plot average lateral velocity
plot([1,vntemp],[minoraxismean_bubbles,minoraxismean_bubbles],'--r') 
xlabel('n_{th} bubble')
ylabel('axis length/mm')
title('major and minor axis of different bubbles')
legend('major axis','minor axis','Location','best')
print('-dtiff',[Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','major and minor axis.tif'])
saveas(gcf,[Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','major and minor axis.fig'])
delete(h)
hold off 

%% plot the equiavlent diameter of different bubbles
clf('reset');close all
h=figure(seed);
 set(h,'visible','on')
plot(xt,diameter_areamean,'k-s','MarkerFaceColor','k','MarkerEdgeColor','k')    
 hold on 
plot(xt,diameter_axismean,'r-o','MarkerFaceColor','r','MarkerEdgeColor','r')
   
% plot each frame velocity
hold on 
plot([1,vntemp],[diameter_areamean_bubbles,diameter_areamean_bubbles],'--k')   % Plot average lateral velocity
hold on
plot([1,vntemp],[diameter_axismean_bubbles,diameter_axismean_bubbles],'--r') 
ylabel('mean diameter of different bubbles/mm')
xlabel('n_{th} bubble')
title('equivalent diameter of different bubbles')
legend('equivalent area diameter','equivalent volume diameter','Location','best')%ellipsoid
print('-dtiff',[Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','equivalent diameter of different bubbles.tif'])
saveas(gcf,[Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','equivalent diameter of different bubbles.fig'])
hold off 
delete(h)
seed=seed+1;
%% plot equivalent area diameter and equivalent volume diameter
    clf('reset');close all
    h=figure(seed);
     set(h,'visible','on');
    scatter(diameter_areamean,diameter_axismean,'k','filled')
    hold on
    diametermeanmax=max([max(diameter_areamean),max(diameter_axismean)]);
    diametermeanmin=min([min(diameter_areamean),min(diameter_axismean)]);
    plot([diametermeanmin,diametermeanmax],[diametermeanmin,diametermeanmax],'--r')   
    xlabel('mean equivalent area diameter/mm')
    ylabel('mean equivalent volume diameter/mm ')
    title('comparison of equivalent mean diameter')
    hold off
    print('-dtiff',[Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','comparison of mean equiavelent area diameter and mean equivalent volume diameter for different bubbles.tif'])
    saveas(gcf,[Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','comparison of mean equiavelent area diameter and mean equivalent volume diameter for different bubbles.fig'])
    delete(h)
    seed=seed+1; 
%% histogram velocity parameters of different bubbles
clf('reset');close all
h=figure(seed);
 set(h,'visible','on');
histogram(vxmean)
title('v_x distribution for different bubbles')
ylabel('frequency')
xlabel('velocity/mm·s^{-1}')
print('-dtiff',[Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','velocity_x distribution of different bubbles.tif'])
saveas(gcf,[Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','velocity_x distribution of different bubbles.fig'])
hold off
seed=seed+1;
clf('reset');close all
h=figure(seed);
 set(h,'visible','on');
histogram(vzmean)
title('v_z distribution for different bubbles')
ylabel('frequency')
xlabel('velocity/mm·s^{-1}')
print('-dtiff',[Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','velocity_z distribution of different bubbles.tif'])
saveas(gcf,[Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','velocity_z distribution of different bubbles.fig'])
hold off
seed=seed+1;
clf('reset');close all
h=figure(seed);
 set(h,'visible','on');
histogram(vmmean)
title('v_m distribution for different bubbles')
ylabel('frequency')
xlabel('velocity/mm·s^{-1}')
print('-dtiff',[Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','velocity_m distribution of different bubbles.tif'])
saveas(gcf,[Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','velocity_m distribution of different bubbles.fig'])
hold off
seed=seed+1;
delete(h)
%% plot diameter distribution
h=figure(seed);
 set(h,'visible','on');
histogram(diameter_areamean)
hold on 
histogram(diameter_axismean)
title('equivalent diameter for different bubbles')
ylabel('frequency')
xlabel('equivalent dianmeter/mm')
legend('equivalent area diameter','equivalent volume diameter','Location','best');
print('-dtiff',[Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','equivalent diameter distribution of different bubbles.tif'])
saveas(gcf,[Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','equivalent diameter distribution of different bubbles.fig'])
hold off
seed=seed+1;
delete(h)

 %% comparison between velocity_x and velocity_z for different bubbles
    clf('reset');close all
    h=figure(seed);
     set(h,'visible','on');
    scatter(vxmean,vzmean,'k','filled')
%     hold on
%     vmeanmax=max([max(vxmean),max(vzmean)]);
%     vmeanmin=min([min(vxmean),min(vzmean)]);
%     plot([vmeanmin,vmeanmax],[vmeanmin,vmeanmax],'--r')   
    xlabel('mean v_x of different bubbles/mm·s^{-1}')
    ylabel('mean v_z of different bubbles/mm·s^{-1}')
    title('comparison of mean velocity_x and mean velocity_z for different bubbles')
    hold off
    print('-dtiff',[Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','comparison of mean velocity_x and mean velocity_z for different bubbles.tif'])
    saveas(gcf,[Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','comparison of mean velocity_x and mean velocity_z for different bubbles.fig'])
    delete(h)
    seed=seed+1;  
 %% comparison of mean minor axis and mean major axis for different bubbles
    clf('reset');close all
    h=figure(seed);
     set(h,'visible','on');
    scatter(minoraxismean,majoraxismean,'k','filled')
    xlabel('mean minor axis/mm')
    ylabel('mean major axis/mm')
    title('comparison of mean minor axis and mean major axis')
    hold off
    print('-dtiff',[Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','comparison of mean minor axis and mean major axis for different bubbles.tif'])
    saveas(gcf,[Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','comparison of mean minor axis and mean major axis for different bubbles.fig'])
    delete(h)
    seed=seed+1;  
%% store the data sheet to files name postprocess*
save([Directorytemp,'PostProcess',num2str(vntemp),'\','FPS',num2str(FPStemp),' ','PostProcess',num2str(vntemp),'.mat'])
 end
end








