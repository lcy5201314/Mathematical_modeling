%% Hopfield�������������䡪������ʶ��
% 
% 
% <html>
% <table border="0" width="600px" id="table1">	<tr>		<td><b><font size="2">�ð�������������</font></b></td>	</tr>	<tr>		<td><span class="comment"><font size="2">1�����˳���פ���ڴ�<a target="_blank" href="http://www.ilovematlab.cn/forum-158-1.html"><font color="#0000FF">���</font></a>���<a target="_blank" href="http://www.ilovematlab.cn/thread-48362-1-1.html"><font color="#0000FF">�ð���</font></a>���ʣ��������ʱش�</font></span></td></tr><tr>	<td><span class="comment"><font size="2">2���˰��������׵Ľ�ѧ��Ƶ�����׵�����������Matlab����</font></span></td>	</tr>	<tr>		<td><span class="comment"><font size="2">		3����������Ϊ�ð����Ĳ������ݣ�Լռ�ð����������ݵ�1/10����</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		4���˰���Ϊԭ��������ת����ע��������<a target="_blank" href="http://www.ilovematlab.cn/">Matlab������̳</a>��<a target="_blank" href="http://www.ilovematlab.cn/forum-158-1.html">��Matlab������30������������</a>����</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		5�����˰��������������о��й��������ǻ�ӭ���������Ҫ��ȣ����ǿ��Ǻ���Լ��ڰ����</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		6������������������Ϊ���壬�鼮��ʵ�����ݿ�����������룬���鼮ʵ�ʷ�������Ϊ׼��</font></span></td>	</tr><tr>		<td><span class="comment"><font size="2">		7�����������������⡢Ԥ����ʽ�ȣ�<a target="_blank" href="http://www.ilovematlab.cn/thread-47939-1-1.html">��������</a>��</font></span></td>	</tr></table>
% </html>
% 
web browser http://www.ilovematlab.cn/thread-60165-1-1.html
%% ��ջ�������
clc
clear
%% ���ݵ���
load data1 array_one
load data2 array_two
%% ѵ��������Ŀ��������
 T=[array_one;array_two]';
%% ��������
 net=newhop(T);
%% ����1��2�Ĵ��������ֵ��󣨹̶�����
load data1_noisy noisy_array_one
load data2_noisy noisy_array_two
%% ����1��2�Ĵ��������ֵ����������
% noisy_array_one=array_one;
% noisy_array_two=array_two;
% for i=1:100
%     a=rand;
%     if a<0.3
%        noisy_array_one(i)=-array_one(i);
%        noisy_array_two(i)=-array_two(i);
%     end
% end
%% ����ʶ��
% identify_one=sim(net,10,[],noisy_array_one');
noisy_one={(noisy_array_one)'};
identify_one=sim(net,{10,10},{},noisy_one);
identify_one{10}';
noisy_two={(noisy_array_two)'};
identify_two=sim(net,{10,10},{},noisy_two);
identify_two{10}';
%% �����ʾ
Array_one=imresize(array_one,20);
subplot(3,2,1)
imshow(Array_one)
title('��׼(����1)') 
Array_two=imresize(array_two,20);
subplot(3,2,2)
imshow(Array_two)
title('��׼(����2)') 
subplot(3,2,3)
Noisy_array_one=imresize(noisy_array_one,20);
imshow(Noisy_array_one)
title('����(����1)') 
subplot(3,2,4)
Noisy_array_two=imresize(noisy_array_two,20);
imshow(Noisy_array_two)
title('����(����2)')
subplot(3,2,5)
imshow(imresize(identify_one{10}',20))
title('ʶ��(����1)')
subplot(3,2,6)
imshow(imresize(identify_two{10}',20))
title('ʶ��(����2)')
web browser http://www.ilovematlab.cn/thread-60165-1-1.html
%%
% 
% <html>
% <table align="center" >	<tr>		<td align="center"><font size="2">��Ȩ���У�</font><a
% href="http://www.ilovematlab.cn/">Matlab������̳</a>&nbsp;&nbsp; <script
% src="http://s3.cnzz.com/stat.php?id=971931&web_id=971931&show=pic" language="JavaScript" ></script>&nbsp;</td>	</tr></table>
% </html>
% 

