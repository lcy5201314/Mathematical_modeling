%% LVQ神经网络的分类——乳腺肿瘤诊断
%
% <html>
% <table border="0" width="600px" id="table1">	<tr>		<td><b><font size="2">该案例作者申明：</font></b></td>	</tr>	<tr>		<td><span class="comment"><font size="2">1：本人长期驻扎在此<a target="_blank" href="http://www.ilovematlab.cn/forum-158-1.html"><font color="#0000FF">板块</font></a>里，对<a target="_blank" href="http://www.ilovematlab.cn/thread-49221-1-1.html"><font color="#0000FF">该案例</font></a>提问，做到有问必答。</font></span></td></tr><tr>	<td><span class="comment"><font size="2">2：此案例有配套的教学视频，配套的完整可运行Matlab程序。</font></span></td>	</tr>	<tr>		<td><span class="comment"><font size="2">		3：以下内容为该案例的部分内容（约占该案例完整内容的1/10）。</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		4：此案例为原创案例，转载请注明出处（<a target="_blank" href="http://www.ilovematlab.cn/">Matlab中文论坛</a>，<a target="_blank" href="http://www.ilovematlab.cn/forum-158-1.html">《Matlab神经网络30个案例分析》</a>）。</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		5：若此案例碰巧与您的研究有关联，我们欢迎您提意见，要求等，我们考虑后可以加在案例里。</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		6：您看到的以下内容为初稿，书籍的实际内容可能有少许出入，以书籍实际发行内容为准。</font></span></td>	</tr><tr>		<td><span class="comment"><font size="2">		7：此书其他常见问题、预定方式等，<a target="_blank" href="http://www.ilovematlab.cn/thread-47939-1-1.html">请点击这里</a>。</font></span></td>	</tr></table>
% </html>
%
web browser http://www.ilovematlab.cn/thread-61656-1-1.html
%% 清空环境变量
clear all
clc
warning off
%% 导入数据
load data.mat
a=randperm(569);
Train=data(a(1:500),:);
Test=data(a(501:end),:);
% 训练数据
P_train=Train(:,3:end)';
Tc_train=Train(:,2)';
T_train=ind2vec(Tc_train);
% 测试数据
P_test=Test(:,3:end)';
Tc_test=Test(:,2)';
%% 创建网络
count_B=length(find(Tc_train==1));
count_M=length(find(Tc_train==2));
rate_B=count_B/500;
rate_M=count_M/500;
net=newlvq(minmax(P_train),20,[rate_B rate_M],0.01,'learnlv1');
% 设置网络参数
net.trainParam.epochs=1000;
net.trainParam.show=10;
net.trainParam.lr=0.1;
net.trainParam.goal=0.1;
%% 训练网络
net=train(net,P_train,T_train);
%% 仿真测试
T_sim=sim(net,P_test);
Tc_sim=vec2ind(T_sim);
result=[Tc_sim;Tc_test]
%% 结果显示
total_B=length(find(data(:,2)==1));
total_M=length(find(data(:,2)==2));
number_B=length(find(Tc_test==1));
number_M=length(find(Tc_test==2));
number_B_sim=length(find(Tc_sim==1 & Tc_test==1));
number_M_sim=length(find(Tc_sim==2 &Tc_test==2));
disp(['病例总数：' num2str(569)...
      '  良性：' num2str(total_B)...
      '  恶性：' num2str(total_M)]);
disp(['训练集病例总数：' num2str(500)...
      '  良性：' num2str(count_B)...
      '  恶性：' num2str(count_M)]);
disp(['测试集病例总数：' num2str(69)...
      '  良性：' num2str(number_B)...
      '  恶性：' num2str(number_M)]);
disp(['良性乳腺肿瘤确诊：' num2str(number_B_sim)...
      '  误诊：' num2str(number_B-number_B_sim)...
      '  确诊率p1=' num2str(number_B_sim/number_B*100) '%']);
disp(['恶性乳腺肿瘤确诊：' num2str(number_M_sim)...
      '  误诊：' num2str(number_M-number_M_sim)...
      '  确诊率p2=' num2str(number_M_sim/number_M*100) '%']);
web browser http://www.ilovematlab.cn/thread-61656-1-1.html  
%%
% 
% <html>
% <table align="center" >	<tr>		<td align="center"><font size="2">版权所有：</font><a
% href="http://www.ilovematlab.cn/">Matlab中文论坛</a>&nbsp;&nbsp; <script
% src="http://s3.cnzz.com/stat.php?id=971931&web_id=971931&show=pic" language="JavaScript" ></script>&nbsp;</td>	</tr></table>
% </html>
% 
