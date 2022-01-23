function chrom=Qgate(chrom,fitness,best,binary)
%% ������ת�ŵ�������
% ����  chrom������ǰ�����ӱ��ر���
%     fitness����Ӧ��ֵ
%        best����ǰ��Ⱥ�����Ÿ���
%      binary�������Ʊ���
% ���  chrom�����º�����ӱ��ر���
sizepop=size(chrom,1)/2;
lenchrom=size(binary,2);
for i=1:sizepop
    for j=1:lenchrom
        A=chrom(2*i-1,j);   % ��
        B=chrom(2*i,j);     % ��
        x=binary(i,j);
        b=best.binary(j);
        if ((x==0)&(b==0))||((x==1)&(b==1))
            delta=0;                  % deltaΪ��ת�ǵĴ�С
            s=0;                        % sΪ��ת�ǵķ��ţ�����ת����
        elseif (x==0)&(b==1)&(fitness(i)<best.fitness)
            delta=0.01*pi;
            if A*B>0
                s=1;
            elseif A*B<0
                s=-1;
            elseif A==0
                s=0;
            elseif B==0
                s=sign(randn);
            end
        elseif (x==0)&(b==1)&(fitness(i)>=best.fitness)
            delta=0.01*pi;
            if A*B>0
                s=-1;
            elseif A*B<0
                s=1;
            elseif A==0
                s=sign(randn);
            elseif B==0
                s=0;
            end
        elseif (x==1)&(b==0)&(fitness(i)<best.fitness)
            delta=0.01*pi;
            if A*B>0
                s=-1;
            elseif A*B<0
                s=1;
            elseif A==0
                s=sign(randn);
            elseif B==0
                s=0;
            end
        elseif (x==1)&(b==0)&(fitness(i)>=best.fitness)
            delta=0.01*pi;
            if A*B>0
                s=1;
            elseif A*B<0
                s=-1;
            elseif A==0
                s=0;
            elseif B==0
                s=sign(randn);
            end
        end
        e=s*delta;       % eΪ��ת��
        U=[cos(e) -sin(e);sin(e) cos(e)];      % ������ת��
        y=U*[A B]';        % yΪ���º������λ
        chrom(2*i-1,j)=y(1);
        chrom(2*i,j)=y(2);
    end
end
