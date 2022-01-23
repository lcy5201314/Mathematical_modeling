%% 输入
% x：一个个体的初始权值和阈值
% P：训练样本输入
% T：训练样本输出
% hiddennum：隐含层神经元数
% P_test:测试样本输入
% T_test:测试样本期望输出
%% 输出
% err：预测样本的预测误差的范数
function [err,T_sim]=BpFunction(x,P,T,hiddennum,P_test,T_test)
inputnum=size(P,1);                             % 输入层神经元个数
% hiddennum=2*inputnum+1;                           % 隐含层神经元个数
outputnum=size(T,1);                                % 输出层神经元个数
%% 数据归一化
[p_train,ps_train]=mapminmax(P,0,1);
p_test=mapminmax('apply',P_test,ps_train);
[t_train,ps_output]=mapminmax(T,0,1);
%% 开始构建BP网络
net=newff(p_train,t_train,hiddennum);               %隐含层为hiddennum个神经元
%设定参数网络参数
net.trainParam.epochs=1000;
net.trainParam.goal=1e-3;
net.trainParam.lr=0.01;
net.trainParam.showwindow=false;                    %高版MATLAB使用 不显示图形框
%% BP神经网络初始权值和阈值
w1num=inputnum*hiddennum;                                           %输入层到隐层的权值个数
w2num=outputnum*hiddennum;                                          %隐含层到输出层的权值个数
% x=2*rand(1,w1num+hiddennum+w2num+outputnum)-1;                      %随即生成权值
W1=x(1:w1num);                                                      %初始输入层到隐含层的权值
B1=x(w1num+1:w1num+hiddennum);                                      %隐层神经元阈值
W2=x(w1num+hiddennum+1:w1num+hiddennum+w2num);                      %隐含层到输出层的权值
B2=x(w1num+hiddennum+w2num+1:w1num+hiddennum+w2num+outputnum);      %输出层阈值
net.iw{1,1}=reshape(W1,hiddennum,inputnum);                         %为神经网络的输入层到隐含层权值赋值
net.lw{2,1}=reshape(W2,outputnum,hiddennum);                        %为神经网络的隐含层到输出层权值赋值
net.b{1}=reshape(B1,hiddennum,1);                                   %为神经网络的隐层神经元阈值赋值
net.b{2}=reshape(B2,outputnum,1);                                   %为神经网络的输出层阈值赋值
%% 开始训练
net = train(net,p_train,t_train);
%% 测试网络
t_sim = sim(net,p_test);
T_sim = mapminmax('reverse',t_sim,ps_output);   %反归一化
err=norm(T_sim-T_test);                         %预测结果与测试结果差的范数，范数越小说明预测得越准确，如果范数为0，说明预测得完全准确

