%% ��ɢHopfield�ķ��ࡪ����У������������
%
% <html>
% <table border="0" width="600px" id="table1">	<tr>		<td><b><font size="2">�ð�������������</font></b></td>	</tr>	<tr>		<td><span class="comment"><font size="2">1�����˳���פ���ڴ�<a target="_blank" href="http://www.ilovematlab.cn/forum-158-1.html"><font color="#0000FF">���</font></a>���<a target="_blank" href="http://www.ilovematlab.cn/thread-49221-1-1.html"><font color="#0000FF">�ð���</font></a>���ʣ��������ʱش�</font></span></td></tr><tr>	<td><span class="comment"><font size="2">2���˰��������׵Ľ�ѧ��Ƶ�����׵�����������Matlab����</font></span></td>	</tr>	<tr>		<td><span class="comment"><font size="2">		3����������Ϊ�ð����Ĳ������ݣ�Լռ�ð����������ݵ�1/10����</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		4���˰���Ϊԭ��������ת����ע��������<a target="_blank" href="http://www.ilovematlab.cn/">Matlab������̳</a>��<a target="_blank" href="http://www.ilovematlab.cn/forum-158-1.html">��Matlab������30������������</a>����</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		5�����˰��������������о��й��������ǻ�ӭ���������Ҫ��ȣ����ǿ��Ǻ���Լ��ڰ����</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		6������������������Ϊ���壬�鼮��ʵ�����ݿ�����������룬���鼮ʵ�ʷ�������Ϊ׼��</font></span></td>	</tr><tr>		<td><span class="comment"><font size="2">		7�����������������⡢Ԥ����ʽ�ȣ�<a target="_blank" href="http://www.ilovematlab.cn/thread-47939-1-1.html">��������</a>��</font></span></td>	</tr></table>
% </html>
%
web browser http://www.ilovematlab.cn/thread-60676-1-1.html
%% ��ջ�������
clear all
clc
%% ��������
load class.mat
%% Ŀ������
T=[class_1 class_2 class_3 class_4 class_5];
%% ��������
net=newhop(T);
%% �������������
load sim.mat
A={[sim_1 sim_2 sim_3 sim_4 sim_5]};
%% �������
Y=sim(net,{25 20},{},A);
%% �����ʾ
Y1=Y{20}(:,1:5)
Y2=Y{20}(:,6:10)
Y3=Y{20}(:,11:15)
Y4=Y{20}(:,16:20)
Y5=Y{20}(:,21:25)
%% ��ͼ
result={T;A{1};Y{20}};
figure
for p=1:3
    for k=1:5 
        subplot(3,5,(p-1)*5+k)
        temp=result{p}(:,(k-1)*5+1:k*5);
        [m,n]=size(temp);
        for i=1:m
            for j=1:n
                if temp(i,j)>0
                   plot(j,m-i,'ko','MarkerFaceColor','k');
                else
                   plot(j,m-i,'ko');
                end
                hold on
            end
        end
        axis([0 6 0 12])
        axis off
        if p==1
           title(['class' num2str(k)])
        elseif p==2
           title(['pre-sim' num2str(k)])
        else
           title(['sim' num2str(k)])
        end
    end                
end
% 
noisy=[1 -1 -1 -1 -1;-1 -1 -1 1 -1;
       -1 1 -1 -1 -1;-1 1 -1 -1 -1;
       1 -1 -1 -1 -1;-1 -1 1 -1 -1;
       -1 -1 -1 1 -1;-1 -1 -1 -1 1;
       -1 1 -1 -1 -1;-1 -1 -1 1 -1;
       -1 -1 1 -1 -1];
y=sim(net,{5 100},{},{noisy});
a=y{100}
web browser http://www.ilovematlab.cn/thread-60676-1-1.html
%%
% 
% <html>
% <table align="center" >	<tr>		<td align="center"><font size="2">��Ȩ���У�</font><a
% href="http://www.ilovematlab.cn/">Matlab������̳</a>&nbsp;&nbsp; <script
% src="http://s3.cnzz.com/stat.php?id=971931&web_id=971931&show=pic" language="JavaScript" ></script>&nbsp;</td>	</tr></table>
% </html>
% 




