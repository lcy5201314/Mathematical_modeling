fs=10^7;
ts=1/fs;%采样时间
fd=0;%doppler
tau=[0,0.2,0.5,1.6,2.3,5.0]*ts*10;%各径时延TU6
pdb=[-3,0,-2,-6,-8,-10]; %各径的平均功率
%tau，pdb应该一一对应，所以位置对应变话，不影响效果？
%tau=[0.2,0,0.5,1.6,2.3,5.0]*ts*10;%实际上对这样导致PathGains变化了，多次平均呢？
%pdb=[0,-3,-2,-6,-8,-10]; 

chan = rayleighchan(ts,fd,tau,pdb);
channel.ResetBeforeFiltering=1;
channel.NormalizePathGains=1;
sigR = filter(sigT,channel);

