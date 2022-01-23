%% ����Hopfield��������Ż��������������Ż�����
%
% <html>
% <table border="0" width="600px" id="table1">	<tr>		<td><b><font size="2">�ð�������������</font></b></td>	</tr>	<tr>		<td><span class="comment"><font size="2">1�����˳���פ���ڴ�<a target="_blank" href="http://www.ilovematlab.cn/forum-158-1.html"><font color="#0000FF">���</font></a>���<a target="_blank" href="http://www.ilovematlab.cn/thread-49221-1-1.html"><font color="#0000FF">�ð���</font></a>���ʣ��������ʱش�</font></span></td></tr><tr>	<td><span class="comment"><font size="2">2���˰��������׵Ľ�ѧ��Ƶ�����׵�����������Matlab����</font></span></td>	</tr>	<tr>		<td><span class="comment"><font size="2">		3����������Ϊ�ð����Ĳ������ݣ�Լռ�ð����������ݵ�1/10����</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		4���˰���Ϊԭ��������ת����ע��������<a target="_blank" href="http://www.ilovematlab.cn/">Matlab������̳</a>��<a target="_blank" href="http://www.ilovematlab.cn/forum-158-1.html">��Matlab������30������������</a>����</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		5�����˰��������������о��й��������ǻ�ӭ���������Ҫ��ȣ����ǿ��Ǻ���Լ��ڰ����</font></span></td>	</tr>		<tr>		<td><span class="comment"><font size="2">		6������������������Ϊ���壬�鼮��ʵ�����ݿ�����������룬���鼮ʵ�ʷ�������Ϊ׼��</font></span></td>	</tr><tr>		<td><span class="comment"><font size="2">		7�����������������⡢Ԥ����ʽ�ȣ�<a target="_blank" href="http://www.ilovematlab.cn/thread-47939-1-1.html">��������</a>��</font></span></td>	</tr></table>
% </html>
%
web browser http://www.ilovematlab.cn/thread-61023-1-1.html
% function main
%% ��ջ�������������ȫ�ֱ���
clear all
clc
global A D
%% �������λ��
load city_location
%% �����໥���м����
distance=dist(citys,citys');
%% ��ʼ������
N=size(citys,1);
A=200;
D=100;
U0=0.1;
step=0.0001;
delta=2*rand(N,N)-1;
U=U0*log(N-1)+delta;
V=(1+tansig(U/U0))/2;
iter_num=10000;
E=zeros(1,iter_num);
%% Ѱ�ŵ���
for k=1:iter_num  
    % ��̬���̼���
    dU=diff_u(V,distance);
    % ������Ԫ״̬����
    U=U+dU*step;
    % �����Ԫ״̬����
    V=(1+tansig(U/U0))/2;
    % ������������
    e=energy(V,distance);
    E(k)=e;  
end
 %% �ж�·����Ч��
[rows,cols]=size(V);
V1=zeros(rows,cols);
[V_max,V_ind]=max(V);
for j=1:cols
    V1(V_ind(j),j)=1;
end
C=sum(V1,1);
R=sum(V1,2);
flag=isequal(C,ones(1,N)) & isequal(R',ones(1,N));
%% �����ʾ
if flag==1
   % �����ʼ·������
   sort_rand=randperm(N);
   citys_rand=citys(sort_rand,:);
   Length_init=dist(citys_rand(1,:),citys_rand(end,:)');
   for i=2:size(citys_rand,1)
       Length_init=Length_init+dist(citys_rand(i-1,:),citys_rand(i,:)');
   end
   % ���Ƴ�ʼ·��
   figure(1)
   plot([citys_rand(:,1);citys_rand(1,1)],[citys_rand(:,2);citys_rand(1,2)],'o-')
   for i=1:length(citys)
       text(citys(i,1),citys(i,2),['   ' num2str(i)])
   end
   text(citys_rand(1,1),citys_rand(1,2),['       ���' ])
   text(citys_rand(end,1),citys_rand(end,2),['       �յ�' ])
   title(['�Ż�ǰ·��(���ȣ�' num2str(Length_init) ')'])
   axis([0 1 0 1])
   grid on
   xlabel('����λ�ú�����')
   ylabel('����λ��������')
   % ��������·������
   [V1_max,V1_ind]=max(V1);
   citys_end=citys(V1_ind,:);
   Length_end=dist(citys_end(1,:),citys_end(end,:)');
   for i=2:size(citys_end,1)
       Length_end=Length_end+dist(citys_end(i-1,:),citys_end(i,:)');
   end
   disp('����·������');V1
   % ��������·��
   figure(2)
   plot([citys_end(:,1);citys_end(1,1)],...
       [citys_end(:,2);citys_end(1,2)],'o-')
   for i=1:length(citys)
       text(citys(i,1),citys(i,2),['  ' num2str(i)])
   end
   text(citys_end(1,1),citys_end(1,2),['       ���' ])
   text(citys_end(end,1),citys_end(end,2),['       �յ�' ])
   title(['�Ż���·��(���ȣ�' num2str(Length_end) ')'])
   axis([0 1 0 1])
   grid on
   xlabel('����λ�ú�����')
   ylabel('����λ��������')
   % �������������仯����
   figure(3)
   plot(1:iter_num,E);
   ylim([0 2000])
   title(['���������仯����(����������' num2str(E(end)) ')']);
   xlabel('��������');
   ylabel('��������');
else
   disp('Ѱ��·����Ч');
end

% %===========================================
% function du=diff_u(V,d)
% global A D
% n=size(V,1);
% sum_x=repmat(sum(V,2)-1,1,n);
% sum_i=repmat(sum(V,1)-1,n,1);
% V_temp=V(:,2:n);
% V_temp=[V_temp V(:,1)];
% sum_d=d*V_temp;
% du=-A*sum_x-A*sum_i-D*sum_d;
% %==========================================
% function E=energy(V,d)
% global A D
% n=size(V,1);
% sum_x=sumsqr(sum(V,2)-1);
% sum_i=sumsqr(sum(V,1)-1);
% V_temp=V(:,2:n);
% V_temp=[V_temp V(:,1)];
% sum_d=d*V_temp;
% sum_d=sum(sum(V.*sum_d));
% E=0.5*(A*sum_x+A*sum_i+D*sum_d);
web browser http://www.ilovematlab.cn/thread-61023-1-1.html
%%
% 
% <html>
% <table align="center" >	<tr>		<td align="center"><font size="2">��Ȩ���У�</font><a
% href="http://www.ilovematlab.cn/">Matlab������̳</a>&nbsp;&nbsp; <script
% src="http://s3.cnzz.com/stat.php?id=971931&web_id=971931&show=pic" language="JavaScript" ></script>&nbsp;</td>	</tr></table>
% </html>
% 
