%% ����8��GRNN������Ԥ�⡪���ڹ���ع�������Ļ�����Ԥ��
%
%
% <html>
% <table border="0" width="600px" id="table1">	<tr>		<td><b><font size="2">�ð�������������</font></b></td>	</tr>	<tr>		<td><span class="comment"><font size="2">1�����˳���פ���ڴ�<a target="_blank" href="http://www.ilovematlab.cn/forum-158-1.html"><font color="#0000FF">���</font></a>���<a target="_blank" href="http://www.ilovematlab.cn/thread-48362-1-1.html"><font color="#0000FF">�ð���</font></a>���ʣ��������ʱش�</font></span></td></tr><tr>	<td><span class="comment"><font size="2">2���˰��������׵Ľ�ѧ��Ƶ�����׵�����������Matlab����</font></span></td>	</tr>	<tr>		<td><span class="comment"><font size="2">		3����������Ϊ�ð����Ĳ������ݣ�Լռ�ð����������ݵ�1/10����</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		4���˰���Ϊԭ��������ת����ע��������<a target="_blank" href="http://www.ilovematlab.cn/">Matlab������̳</a>��<a target="_blank" href="http://www.ilovematlab.cn/forum-158-1.html">��Matlab������30������������</a>����</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		5�����˰��������������о��й��������ǻ�ӭ���������Ҫ��ȣ����ǿ��Ǻ���Լ��ڰ����</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		6������������������Ϊ���壬�鼮��ʵ�����ݿ�����������룬���鼮ʵ�ʷ�������Ϊ׼��</font></span></td>	</tr><tr>		<td><span class="comment"><font size="2">		7�����������������⡢Ԥ����ʽ�ȣ�<a target="_blank" href="http://www.ilovematlab.cn/thread-47939-1-1.html">��������</a>��</font></span></td>	</tr></table>
% </html>
%
%% ��ջ�������
clc;
clear all
close all
nntwarn off;

%% ��������
load data;
% �������ݲ������ݷֳ�ѵ����Ԥ������
p_train=p(1:12,:);
t_train=t(1:12,:);
p_test=p(13,:);
t_test=t(13,:);
%% ������֤
desired_spread=[];
mse_max=10e20;
desired_input=[];
desired_output=[];
result_perfp=[];
indices = crossvalind('Kfold',length(p_train),4);
h=waitbar(0,'����Ѱ�����Ż�����....')
k=1;
for i = 1:4
    perfp=[];
    disp(['����Ϊ��',num2str(i),'�ν�����֤���'])
    test = (indices == i); train = ~test;
    p_cv_train=p_train(train,:);
    t_cv_train=t_train(train,:);
    p_cv_test=p_train(test,:);
    t_cv_test=t_train(test,:);
    p_cv_train=p_cv_train';
    t_cv_train=t_cv_train';
    p_cv_test= p_cv_test';
    t_cv_test= t_cv_test';
    [p_cv_train,minp,maxp,t_cv_train,mint,maxt]=premnmx(p_cv_train,t_cv_train);
    p_cv_test=tramnmx(p_cv_test,minp,maxp);
    for spread=0.1:0.1:2;
        net=newgrnn(p_cv_train,t_cv_train,spread);
        waitbar(k/80,h);
        disp(['��ǰspreadֵΪ', num2str(spread)]);
        test_Out=sim(net,p_cv_test);
        test_Out=postmnmx(test_Out,mint,maxt);
        error=t_cv_test-test_Out;
        disp(['��ǰ�����mseΪ',num2str(mse(error))])
        perfp=[perfp mse(error)];
        if mse(error)<mse_max
            mse_max=mse(error);
            desired_spread=spread;
            desired_input=p_cv_train;
            desired_output=t_cv_train;
        end
        k=k+1;
    end
    result_perfp(i,:)=perfp;
end;
close(h)
disp(['���spreadֵΪ',num2str(desired_spread)])
disp(['��ʱ�������ֵΪ'])
desired_input
disp(['��ʱ������ֵΪ'])
desired_output
%% ������ѷ�������GRNN����
net=newgrnn(desired_input,desired_output,desired_spread);
p_test=p_test';
p_test=tramnmx(p_test,minp,maxp);
grnn_prediction_result=sim(net,p_test);
grnn_prediction_result=postmnmx(grnn_prediction_result,mint,maxt);
grnn_error=t_test-grnn_prediction_result';
disp(['GRNN��������������Ԥ������Ϊ',num2str(abs(grnn_error))])
save best desired_input desired_output p_test t_test grnn_error mint maxt

web browser http://www.ilovematlab.cn/thread-65171-1-1.html
%%
%
% <html>
% <table align="center" >	<tr>		<td align="center"><font size="2">��Ȩ���У�</font><a
% href="http://www.ilovematlab.cn/">Matlab������̳</a>&nbsp;&nbsp; <script
% src="http://s3.cnzz.com/stat.php?id=971931&web_id=971931&show=pic" language="JavaScript" ></script>&nbsp;</td>	</tr></table>
% </html>
%



     
 
 
 
 
 
 
 