function flag = wj_wfs_con(D)
%�����ڽӾ��������ͨ����1�����򷵻��㡣

L = size(D,1);
visited = [1];
while(length(visited) < L)
    D(visited, visited) = 0;
    pos = D(visited,:);
    if sum(sum(pos)) == 0
        flag = 0;
        return;
    else
        for i = 1:size(pos,1)
            b = find(pos(i,:));
            visited = union(visited, b)
            D(visited, visited) = 0;
        end
    end
end
flag = 1;

