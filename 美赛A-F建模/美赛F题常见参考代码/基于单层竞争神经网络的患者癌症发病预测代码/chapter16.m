%% ����16�����㾺������������ݷ��ࡪ���߰�֢����Ԥ��
% 
% 
% <html>
% <table border="0" width="600px" id="table1">	<tr>		<td><b><font size="2">�ð�������������</font></b></td>	</tr>	<tr>		<td><span class="comment"><font size="2">1�����˳���פ���ڴ�<a target="_blank" href="http://www.ilovematlab.cn/forum-158-1.html"><font color="#0000FF">���</font></a>���<a target="_blank" href="http://www.ilovematlab.cn/thread-48362-1-1.html"><font color="#0000FF">�ð���</font></a>���ʣ��������ʱش�</font></span></td></tr><tr>	<td><span class="comment"><font size="2">2���˰��������׵Ľ�ѧ��Ƶ�����׵�����������Matlab����</font></span></td>	</tr>	<tr>		<td><span class="comment"><font size="2">		3����������Ϊ�ð����Ĳ������ݣ�Լռ�ð����������ݵ�1/10����</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		4���˰���Ϊԭ��������ת����ע��������<a target="_blank" href="http://www.ilovematlab.cn/">Matlab������̳</a>��<a target="_blank" href="http://www.ilovematlab.cn/forum-158-1.html">��Matlab������30������������</a>����</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		5�����˰��������������о��й��������ǻ�ӭ���������Ҫ��ȣ����ǿ��Ǻ���Լ��ڰ����</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		6������������������Ϊ���壬�鼮��ʵ�����ݿ�����������룬���鼮ʵ�ʷ�������Ϊ׼��</font></span></td>	</tr><tr>		<td><span class="comment"><font size="2">		7�����������������⡢Ԥ����ʽ�ȣ�<a target="_blank" href="http://www.ilovematlab.cn/thread-47939-1-1.html">��������</a>��</font></span></td>	</tr></table>
% </html>
% 


%% ��ջ�������
clc
clear

%% ¼����������
% �������ݲ������ݷֳ�ѵ����Ԥ������
load gene.mat;
data=gene;
P=data(1:40,:);
T=data(41:60,:);

% ת�ú����������������ʽ
P=P';
T=T';
% ȡ����Ԫ�ص����ֵ����СֵQ��
Q=minmax(P);

%% ���罨����ѵ��
% ����newc( )������������磺2�����������Ԫ������Ҳ����Ҫ����ĸ�����0.1����ѧϰ���ʡ�
net=newc(Q,2,0.1)

% ��ʼ�����缰�趨���������
net=init(net);
net.trainparam.epochs=20;
% ѵ�����磺
net=train(net,P);


%% �����Ч����֤

% ��ԭ���ݻش�����������Ч����
a=sim(net,P);
ac=vec2ind(a)

% ����ʹ���˱任����vec2ind()�����ڽ���ֵ������任���±�����������õĸ�ʽΪ��
%  ind=vec2ind(vec)
% ���У�
% vec��Ϊm��n�е���������x��x�е�ÿ��������i��������һ��1�⣬����Ԫ�ؾ�Ϊ0��
% ind��Ϊn��Ԫ��ֵΪ1���ڵ����±�ֵ���ɵ�һ����������



%% �����������Ԥ��
% ���潫��20�����ݴ���������ģ���У��۲����������
% sim( )�����������
Y=sim(net,T)
yc=vec2ind(Y)

web browser http://www.ilovematlab.cn/viewthread.php?tid=60656
%%
% 
% <html>
% <table align="center" >	<tr>		<td align="center"><font size="2">��Ȩ���У�</font><a
% href="http://www.ilovematlab.cn/">Matlab������̳</a>&nbsp;&nbsp; <script
% src="http://s3.cnzz.com/stat.php?id=971931&web_id=971931&show=pic" language="JavaScript" ></script>&nbsp;</td>	</tr></table>
% </html>
% 

