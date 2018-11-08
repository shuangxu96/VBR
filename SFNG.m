function SFNet = SFNG(Nodes, mlinks)
%seed =ER(10,0.5);
seed =[0 1 0 0 1;1 0 0 1 0;0 0 0 1 0;0 1 1 0 0;1 0 0 0 0];
if mlinks>5
    seed = blkdiag(seed,seed);
elseif mlinks>10
    seed = blkdiag(seed,seed,seed);
elseif mlinks>15
    seed = blkdiag(seed,seed,seed,seed);
end
seed = full(seed);
pos = length(seed);

%if (Nodes < pos) || (mlinks > pos) || (pos < 1) || (size(size(seed)) ~= 2) || (mlinks < 1) || (seed ~= seed') || (sum(diag(seed)) ~= 0)
%    error('invalid parameter value(s)');
%end

%if mlinks > 5 || Nodes > 15000 || pos > 15000
%    warning('Abnormally large value(s) may cause long processing time');
%end

rand('state',sum(100*clock));

Net = zeros(Nodes, Nodes, 'single');
Net(1:pos,1:pos) = seed;
sumlinks = sum(sum(Net));

while pos < Nodes
    pos = pos + 1;
    linkage = 0;
    while linkage ~= mlinks
        rnode = ceil(rand * pos);
        deg = sum(Net(:,rnode)) * 2;
        rlink = rand * 1;
        if rlink < deg / sumlinks && Net(pos,rnode) ~= 1 && Net(rnode,pos) ~= 1
            Net(pos,rnode) = 1;
            Net(rnode,pos) = 1;
            linkage = linkage + 1;
            sumlinks = sumlinks + 2;
        end
    end
end
index = find(diag(Net)==1);
if ~isempty(index)
    Net(index,index) = 0;
end
clear Nodes deg linkage pos rlink rnode sumlinks mlinks
SFNet = Net;