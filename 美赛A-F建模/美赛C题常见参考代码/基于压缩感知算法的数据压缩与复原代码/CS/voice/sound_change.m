%语音信号采集、处理、传输、显示 
clear;
clc;

N=4096; 
M=25600; 
[y,fs,bite]=wavread('wrd.wav',[N,M]);%采集语音信号
sound(y,fs,bite);
%对语音信号加噪声 
s=wavread('wrd.wav',M); 
s=s+0.03*randn(1,M)';%幅度为0.03 
sound(s,fs,bite);
%用小波db3对s进行2层分解 
level=2; 
[c,l]=wavedec(s,level,'db3'); 
%选用全局阈值进行信号增强处理 
thr=2; 
[sd,csd,lsd,perf0,perf12]=wdencmp('gbl',c,l,'db3',level,thr,'h',1); 
% prompt={'Enter matrix size:'}; 
% dlg_title='Input for peaks function'; 
% num_lines=1; 
% def={'1'};%用户界面对话框，选择回放语音信号的类型 
% answer=inputdlg(prompt,dlg_title,num_lines,def); 
%m=answer; 
%if answer==10 
 sound(sd,fs,bite);%回放语音信号y/s/sd 
%else 
 % sound(s,fs,bite); 
%end 
     
%在图像中显示 
subplot(3,2,1);plot(y);title('原音波形') 
subplot(3,2,3);plot(s);title('加噪声后的波形') 
subplot(3,2,5);plot(sd);title('去噪声后的波形') 
Y=fft(y,4096);%对语音信号进行fft 
subplot(3,2,2);plot(abs(Y));title('语音信号的频谱'); 
S=fft(s,4096); 
subplot(3,2,4);plot(abs(S));title('加噪后语音信号的频谱') 
SD=fft(sd,4096); 
subplot(3,2,6);plot(abs(SD));title('去噪后语音信号的频谱') 
%%%%%%%%%%%%%%%%%%%%%%%%%详细说明：录制一段个人自己的语音信号。对录制的信号进行采样；
%画出采样后语音信号的时域波形和频谱图；对语音信号进行加噪和去噪处理，画出滤
%波后信号的时域波形和频谱，并对滤波前后的信号进行对比，分析信号的变化；回放语音信号；实现快录慢放、慢录快放等功能。-Record a person' s own voice signal. Of the recorded signal sampling draw sample after the speech signal time-domain waveform and frequency spectrum of voice signals, noise and de-noising processing, draw Filtered time-domain signal waveform and spectrum, and filtering the signal before and after comparative analysis of signal changes playback voice signal realize quickly recorded slow release, slow release recorded faster functions.