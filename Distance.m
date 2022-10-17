function [D]=Distance(M)
Temp=ones(size(M));
L{1}=M;
S=M;
i=1;
while 1
    Temp=logical(L{i}*M)-S;
    Temp(Temp<0)=0;
    Temp=Temp-diag(diag(Temp));
    Temp=logical(Temp);
    i=i+1;
    L{i}=Temp;
    S=S+i*L{i};
    if ~any(any(L{i}))
        break;
    end
    %disp(num2str(i));
end
D=S;
end