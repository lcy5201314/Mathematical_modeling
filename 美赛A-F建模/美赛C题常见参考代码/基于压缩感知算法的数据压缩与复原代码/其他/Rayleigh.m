fs=10^7;
ts=1/fs;%����ʱ��
fd=0;%doppler
tau=[0,0.2,0.5,1.6,2.3,5.0]*ts*10;%����ʱ��TU6
pdb=[-3,0,-2,-6,-8,-10]; %������ƽ������
%tau��pdbӦ��һһ��Ӧ������λ�ö�Ӧ�仰����Ӱ��Ч����
%tau=[0.2,0,0.5,1.6,2.3,5.0]*ts*10;%ʵ���϶���������PathGains�仯�ˣ����ƽ���أ�
%pdb=[0,-3,-2,-6,-8,-10]; 

chan = rayleighchan(ts,fd,tau,pdb);
channel.ResetBeforeFiltering=1;
channel.NormalizePathGains=1;
sigR = filter(sigT,channel);

