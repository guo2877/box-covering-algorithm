function [Boxes,Nb] = SM30(Net,Lb,n)
% INPUTS: adjacency matrix, box size,sampling density 
% OUTPUTS:box coverage, the number of boxes
%！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
tic;
D=Distance(Net); 
D(D>=Lb)=0;
D(D~=0)=1;

DB=maximal_box_Sampling(D,n);   
freq=sum(DB,1);   
Boxes=zeros(1,size(DB,2));   
Color=1; 

r_sums=sum(DB,2);  
[~, index]=sort(r_sums);  
DB=DB(index,:);    
while ~isequal(freq,zeros(1,size(D,2)))
    cnode=find(DB(1,:)==1);
    fvalue=freq(cnode);  
    if size(find(fvalue==1),2)==0   
        freq(cnode)=freq(cnode)-1;  
        DB(1,:)=[];
    else
        Boxes(cnode)=Color;
        Color=Color+1;
        freq(cnode)=0;  
        DB(:,cnode)=0;
        DB(all(DB==0,2),:)=[];  
        r_sums=sum(DB,2);  
        [~, index]=sort(r_sums);   
        DB=DB(index,:);    
    end
end
Nb=max(Boxes);
toc;
end

function [MC] = maximal_box_Sampling(Net,n)        
Scale=size(Net,2);       
S=1:Scale;           
MC_cell=cell(1,1);  
j=0;
while ~isempty(S)
    j=j+1;
    u=S(randperm(length(S),1));   
    [R]=ball_of_seed(Net,u);
    S=setdiff(S,R);
    MC_cell{j}=R;     
end
disp(num2str(j));  
runs=j*n;  
parfor i=j+1:runs   
    u=randperm(Scale,1);  
    [R]=ball_of_seed(Net,u);
    MC_cell{i}=R;    
end
MC=zeros(runs,Scale);
for i=1:runs
    MC(i,MC_cell{i})=1;   
end 
MC=unique(MC,'rows','stable');   
disp(num2str(runs));  
end

function [R]=ball_of_seed(Net,u)
R=u;
P=find(Net(u,:)==1);  
while ~isempty(P)     
    ui=P(randperm(length(P),1));   
    P=setdiff(P,ui);     
    R=union(R,ui);       
    P=intersect(P,find(Net(ui,:)==1));  
end
end