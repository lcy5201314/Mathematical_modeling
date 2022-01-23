%% SVM����������ݷ���Ԥ��----��������Ѿ�����ʶ��
% 
% 
% <html>
% <table border="0" width="600px" id="table1">	<tr>		<td><b><font size="2">�ð�������������</font></b></td>	</tr>	<tr>		<td><span class="comment"><font size="2">1�����˳���פ���ڴ�<a target="_blank" href="http://www.ilovematlab.cn/forum-158-1.html"><font color="#0000FF">���</font></a>���<a target="_blank" href="http://www.ilovematlab.cn/thread-48362-1-1.html"><font color="#0000FF">�ð���</font></a>���ʣ��������ʱش�</font></span></td></tr><tr>	<td><span class="comment"><font size="2">2���˰��������׵Ľ�ѧ��Ƶ�����׵�����������Matlab����</font></span></td>	</tr>	<tr>		<td><span class="comment"><font size="2">		3����������Ϊ�ð����Ĳ������ݣ�Լռ�ð����������ݵ�1/10����</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		4���˰���Ϊԭ��������ת����ע��������<a target="_blank" href="http://www.ilovematlab.cn/">Matlab������̳</a>��<a target="_blank" href="http://www.ilovematlab.cn/forum-158-1.html">��Matlab������30������������</a>����</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		5�����˰��������������о��й��������ǻ�ӭ���������Ҫ��ȣ����ǿ��Ǻ���Լ��ڰ����</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		6������������������Ϊ���壬�鼮��ʵ�����ݿ�����������룬���鼮ʵ�ʷ�������Ϊ׼��</font></span></td>	</tr><tr>		<td><span class="comment"><font size="2">		7�����������������⡢Ԥ����ʽ�ȣ�<a target="_blank" href="http://www.ilovematlab.cn/thread-47939-1-1.html">��������</a>��</font></span></td>	</tr></table>
% </html>
% 


%% ��ջ�������
close all;
clear;
clc;
format compact;
%% ������ȡ

% �����������wine,���а���������Ϊclassnumber = 3,wine:178*13�ľ���,wine_labes:178*1��������
load chapter12_wine.mat;

% �����������ݵ�box���ӻ�ͼ
figure;
boxplot(wine,'orientation','horizontal','labels',categories);
title('wine���ݵ�box���ӻ�ͼ','FontSize',12);
xlabel('����ֵ','FontSize',12);
grid on;

% �����������ݵķ�ά���ӻ�ͼ
figure
subplot(3,5,1);
hold on
for run = 1:178
    plot(run,wine_labels(run),'*');
end
xlabel('����','FontSize',10);
ylabel('����ǩ','FontSize',10);
title('class','FontSize',10);
for run = 2:14
    subplot(3,5,run);
    hold on;
    str = ['attrib ',num2str(run-1)];
    for i = 1:178
        plot(i,wine(i,run-1),'*');
    end
    xlabel('����','FontSize',10);
    ylabel('����ֵ','FontSize',10);
    title(str,'FontSize',10);
end

% ѡ��ѵ�����Ͳ��Լ�

% ����һ���1-30,�ڶ����60-95,�������131-153��Ϊѵ����
train_wine = [wine(1:30,:);wine(60:95,:);wine(131:153,:)];
% ��Ӧ��ѵ�����ı�ǩҲҪ�������
train_wine_labels = [wine_labels(1:30);wine_labels(60:95);wine_labels(131:153)];
% ����һ���31-59,�ڶ����96-130,�������154-178��Ϊ���Լ�
test_wine = [wine(31:59,:);wine(96:130,:);wine(154:178,:)];
% ��Ӧ�Ĳ��Լ��ı�ǩҲҪ�������
test_wine_labels = [wine_labels(31:59);wine_labels(96:130);wine_labels(154:178)];

%% ����Ԥ����
% ����Ԥ����,��ѵ�����Ͳ��Լ���һ����[0,1]����

[mtrain,ntrain] = size(train_wine);
[mtest,ntest] = size(test_wine);

dataset = [train_wine;test_wine];
% mapminmaxΪMATLAB�Դ��Ĺ�һ������
[dataset_scale,ps] = mapminmax(dataset',0,1);
dataset_scale = dataset_scale';

train_wine = dataset_scale(1:mtrain,:);
test_wine = dataset_scale( (mtrain+1):(mtrain+mtest),: );
%% SVM����ѵ��
model = svmtrain(train_wine_labels, train_wine, '-c 2 -g 1');

%% SVM����Ԥ��
[predict_label, accuracy] = svmpredict(test_wine_labels, test_wine, model);

%% �������

% ���Լ���ʵ�ʷ����Ԥ�����ͼ
% ͨ��ͼ���Կ���ֻ��һ�����������Ǳ���ֵ�
figure;
hold on;
plot(test_wine_labels,'o');
plot(predict_label,'r*');
xlabel('���Լ�����','FontSize',12);
ylabel('����ǩ','FontSize',12);
legend('ʵ�ʲ��Լ�����','Ԥ����Լ�����');
title('���Լ���ʵ�ʷ����Ԥ�����ͼ','FontSize',12);
grid on;

%%
% 
% <html>
% <table align="center" >	<tr>		<td align="center"><font size="2">��Ȩ���У�</font><a
% href="http://www.ilovematlab.cn/">Matlab������̳</a>&nbsp;&nbsp; <script
% src="http://s3.cnzz.com/stat.php?id=971931&web_id=971931&show=pic" language="JavaScript" ></script>&nbsp;</td>	</tr></table>
% </html>
% 

