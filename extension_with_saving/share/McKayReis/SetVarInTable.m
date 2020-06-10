function  SetVarInTable( beta, b , B, R, h,  lambda, u, wage, MU, stdlogy, stdu, welfareSS)

load('results/WorkTable.mat')
WorkTableIndices;


row = (abs(WorkTable(:,ibeta) - beta) < 1e-6) & (abs(WorkTable(:,ib) - b) < 1e-6);
if sum(row) == 0
    row = size(WorkTable,1)+1;
    WorkTable = [WorkTable; zeros(1,size(WorkTable,2))];
    WorkTable(row,1:2) = [beta b];
else
    assert(sum(row) == 1);
end


rowvar = WorkTable(row,:);


if exist('B','var') && ~isempty(B)
    rowvar(iB)=B;
end
if exist('R','var') && ~isempty(R)
    rowvar(iR) = R;
end
if exist('h','var') && ~isempty(h)
    rowvar(ih) = h;
end
if exist('lambda','var') && ~isempty(lambda)
    rowvar(ilambda) = lambda;
end
if exist('u','var') && ~isempty(u)
    rowvar(iu) = u;
end
if exist('wage','var') && ~isempty(wage)
    rowvar(iwage) = wage;
end
if exist('MU','var') && ~isempty(MU)
    rowvar(iMU) = MU;
end
if exist('stdlogy','var') && ~isempty(stdlogy)
    rowvar(istdlogy) = stdlogy;
end
if exist('stdu','var') && ~isempty(stdu)
    rowvar(istdu) = stdu;
end
if exist('welfareSS','var') && ~isempty(welfareSS)
    rowvar(iwelfareSS) = welfareSS;
end

WorkTable(row,:) = rowvar;

save('results/WorkTable.mat','WorkTable');

end

