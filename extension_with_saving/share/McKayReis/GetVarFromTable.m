function [ B, R, h, q, lambda, u, wage, MU, stdlogy ] = GetVarFromTable( beta, b, Strict )

if ~exist('Strict','var')
    Strict = true;
end

load('results/WorkTable.mat')
WorkTableIndices;

if Strict
    row = (abs(WorkTable(:,ibeta) - beta) < 1e-6) & (abs(WorkTable(:,ib) - b) < 1e-6);
    assert(sum(row) == 1);
else
    [~,row] = min( 1000*abs(WorkTable(:,ibeta) - beta) +  abs(WorkTable(:,ib) - b));
end
 
row = WorkTable(row,:);

B = row(iB);
R = row(iR);
h = row(ih);
q = nan;
lambda = row(ilambda);
u = row(iu);
wage = row(iwage);
MU = row(iMU);
stdlogy = row(istdlogy);

end

