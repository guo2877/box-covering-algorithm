function [Boxes,Nb] = OBCA(Net,Lb)
% INPUTS: adjacency matrix, box size
% OUTPUTS:box coverage, the number of boxes
%！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
tic;
D=Distance(Net); 
D(D>=Lb)=0;
D(D~=0)=1;   
Scale=size(D,1);   
Boxes=[];    
freq =zeros(1,Scale);         
[deg,~,~]=degrees(Net);   
[~, index]=sort(deg);     
runs=1;
S=index;

for n=1:Scale
    seed=index(n);
    if ~ismember(seed,S)
        continue;
    else
        box=ball_of_seed(D,seed,S);
        d=find(box==1);
        for i=1:size(d,2)
            S(S==d(i))=[];
        end
        Boxes(runs,:)=box;  
        runs=runs+1;
        freq(1,box==1)=freq(1,box==1)+1;  
    end
end

disp(num2str(runs)); 
N=size(Boxes,1);
d_box=[];
for n=1:N
    redundancy=0;
    cnode=find(Boxes(n,:)==1);
    num_node=size(cnode,2);
    for i=1:num_node
        if freq(cnode(i))>1
            redundancy=redundancy+1;
        end
    end
    if ~isequal(redundancy,num_node)
        continue;
    else
        d_box=[d_box,n];
        for i=1:num_node
            freq(cnode(i))=freq(cnode(i))-1;
        end
    end
end

Boxes(d_box,:)=[];
Nb=size(Boxes,1);
toc;
end

function [box]=ball_of_seed(D,seed,S)
box=zeros(1,size(D,1));   
box(1,seed)=1;
P=find(D(seed,:)==1);  
while ~isempty(P)
    pre_nodes=intersect(S,P);
    if ~isequal(size(pre_nodes,2),0)
        ui=pre_nodes(randperm(length(pre_nodes),1)); 
    else
        ui=P(randperm(length(P),1)); 
    end
    P=setdiff(P,ui);    
    box(1,ui)=1;       
    P=intersect(P,find(D(ui,:)==1));  
end
end

