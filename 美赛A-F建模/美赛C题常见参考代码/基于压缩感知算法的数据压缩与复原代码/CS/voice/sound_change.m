%�����źŲɼ����������䡢��ʾ 
clear;
clc;

N=4096; 
M=25600; 
[y,fs,bite]=wavread('wrd.wav',[N,M]);%�ɼ������ź�
sound(y,fs,bite);
%�������źż����� 
s=wavread('wrd.wav',M); 
s=s+0.03*randn(1,M)';%����Ϊ0.03 
sound(s,fs,bite);
%��С��db3��s����2��ֽ� 
level=2; 
[c,l]=wavedec(s,level,'db3'); 
%ѡ��ȫ����ֵ�����ź���ǿ���� 
thr=2; 
[sd,csd,lsd,perf0,perf12]=wdencmp('gbl',c,l,'db3',level,thr,'h',1); 
% prompt={'Enter matrix size:'}; 
% dlg_title='Input for peaks function'; 
% num_lines=1; 
% def={'1'};%�û�����Ի���ѡ��ط������źŵ����� 
% answer=inputdlg(prompt,dlg_title,num_lines,def); 
%m=answer; 
%if answer==10 
 sound(sd,fs,bite);%�ط������ź�y/s/sd 
%else 
 % sound(s,fs,bite); 
%end 
     
%��ͼ������ʾ 
subplot(3,2,1);plot(y);title('ԭ������') 
subplot(3,2,3);plot(s);title('��������Ĳ���') 
subplot(3,2,5);plot(sd);title('ȥ������Ĳ���') 
Y=fft(y,4096);%�������źŽ���fft 
subplot(3,2,2);plot(abs(Y));title('�����źŵ�Ƶ��'); 
S=fft(s,4096); 
subplot(3,2,4);plot(abs(S));title('����������źŵ�Ƶ��') 
SD=fft(sd,4096); 
subplot(3,2,6);plot(abs(SD));title('ȥ��������źŵ�Ƶ��') 
%%%%%%%%%%%%%%%%%%%%%%%%%��ϸ˵����¼��һ�θ����Լ��������źš���¼�Ƶ��źŽ��в�����
%���������������źŵ�ʱ���κ�Ƶ��ͼ���������źŽ��м����ȥ�봦��������
%�����źŵ�ʱ���κ�Ƶ�ף������˲�ǰ����źŽ��жԱȣ������źŵı仯���ط������źţ�ʵ�ֿ�¼���š���¼��ŵȹ��ܡ�-Record a person' s own voice signal. Of the recorded signal sampling draw sample after the speech signal time-domain waveform and frequency spectrum of voice signals, noise and de-noising processing, draw Filtered time-domain signal waveform and spectrum, and filtering the signal before and after comparative analysis of signal changes playback voice signal realize quickly recorded slow release, slow release recorded faster functions.