%% Hopfield�������������䡪������ʶ��
% 
% 
% <html>
% <table border="0" width="600px" id="table1">	<tr>		<td><b><font size="2">�ð�������������</font></b></td>	</tr>	<tr>		<td><span class="comment"><font size="2">1�����˳���פ���ڴ�<a target="_blank" href="http://www.ilovematlab.cn/forum-158-1.html"><font color="#0000FF">���</font></a>���<a target="_blank" href="http://www.ilovematlab.cn/thread-48362-1-1.html"><font color="#0000FF">�ð���</font></a>���ʣ��������ʱش�</font></span></td></tr><tr>	<td><span class="comment"><font size="2">2���˰��������׵Ľ�ѧ��Ƶ�����׵�����������Matlab����</font></span></td>	</tr>	<tr>		<td><span class="comment"><font size="2">		3����������Ϊ�ð����Ĳ������ݣ�Լռ�ð����������ݵ�1/10����</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		4���˰���Ϊԭ��������ת����ע��������<a target="_blank" href="http://www.ilovematlab.cn/">Matlab������̳</a>��<a target="_blank" href="http://www.ilovematlab.cn/forum-158-1.html">��Matlab������30������������</a>����</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		5�����˰��������������о��й��������ǻ�ӭ���������Ҫ��ȣ����ǿ��Ǻ���Լ��ڰ����</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		6������������������Ϊ���壬�鼮��ʵ�����ݿ�����������룬���鼮ʵ�ʷ�������Ϊ׼��</font></span></td>	</tr><tr>		<td><span class="comment"><font size="2">		7�����������������⡢Ԥ����ʽ�ȣ�<a target="_blank" href="http://www.ilovematlab.cn/thread-47939-1-1.html">��������</a>��</font></span></td>	</tr></table>
% </html>
% 
web browser http://www.ilovematlab.cn/thread-60165-1-1.html
%% �����������
clear all
clc
%% �������ģʽ
load data1.mat
T=array_one;
%% ���������Ȩϵ������
[m,n]=size(T);
w=zeros(m);
for i=1:n
    w=w+T(:,i)*T(:,i)'-eye(m);
end
%% ���������ģʽ
noisy_array=T;
for i=1:100
    a=rand;
    if a<0
       noisy_array(i)=-T(i);
    end
end
%% ��������
v0=noisy_array;
v=zeros(m,n);
for k=1:5
    for i=1:m
        v(i,:)=sign(w(i,:)*v0);
    end
    v0=v;
end
%% ��ͼ
subplot(3,1,1)
t=imresize(T,20);
imshow(t)
title('��׼')
subplot(3,1,2)
Noisy_array=imresize(noisy_array,20);
imshow(Noisy_array)
title('����')
subplot(3,1,3)
V=imresize(v,20);
imshow(V)
title('ʶ��')
web browser http://www.ilovematlab.cn/thread-60165-1-1.html
%%
% 
% <html>
% <table align="center" >	<tr>		<td align="center"><font size="2">��Ȩ���У�</font><a
% href="http://www.ilovematlab.cn/">Matlab������̳</a>&nbsp;&nbsp; <script
% src="http://s3.cnzz.com/stat.php?id=971931&web_id=971931&show=pic" language="JavaScript" ></script>&nbsp;</td>	</tr></table>
% </html>
% 