function caliinc(npp)
global Params;

Params.npp = npp;

x = [0.9;1.1];
[Params.xY,iflag] = broydn(@inciid_z,x,[1e-8,0,0]);
assert(iflag==0);
[resid,Params.wY] =  inciid_z(Params.xY);


if(Params.npp==1)
  Params.ypp = 1;
  Params.transpp = 1;
elseif(Params.npp==2)
  %Params.ypp = [1 1];
  Params.ypp = [0.5 1.5];
  Params.transpp = [0.9 0.1;0.1 0.9];
elseif(Params.npp>10)


  n = Params.npp;
  assert(mod(n,2)==1);
  nmed = (n+1)/2;  %index of median


  
  Params.ypp = exp(linspace(-1.2,1.2,n));
  % make biggest income group very big:
  % Params.ypp(n) = Params.ypp(n)*(Params.ypp(2)/Params.ypp(1))^(fac-1);

  pSwitch = sqrt(0.015)/diff(log(Params.ypp(1:2)))/Params.freq;
  pSwitch2 = pSwitch/2;
  x = [pSwitch;pSwitch2;1.5];
  [x,iflag] = broydn(@cali_z,x,[1e-7,0,1],nmed);
  [resid,Pi] = cali_z(x,nmed);
  
  % make it global:
  Params.transpp = Pi;
  Params.distpp = invdistr(Pi')';
else
  error('wrong npp');
end




function [resid,Pi] = cali_z(x,nmed)
global Params;

pSwitch = x(1);
pSwitch2 = x(2);
expo = x(3);
if(expo<1)
  resid = 1e100;
  return;
end


% transition probability of survivors
PiHH = make_pi(pSwitch,pSwitch2,Params.npp); 
Pi4 = PiHH^Params.freq;
pdeath = 0.02/Params.freq;
PiDeath = eye(Params.npp)*(1-pdeath);
PiDeath(:,nmed) = PiDeath(:,nmed) + pdeath;
Pi = PiDeath*PiHH;
Params.distpp = invdistr(Pi')';

frac1 = 0.1476;
Params.ypp = exp(linspace(-1.2,1.2,Params.npp));
nbig = ceil(4*Params.npp/5);  %index of median
Params.ypp(nbig:end) = Params.ypp(nbig:end).^expo;
ytop = 100*frac1/Params.xY(end);
ytarget = (1-frac1-Params.distpp(end)*Params.wY(1)*ytop*Params.xY(1)); %/(0.99-Params.distpp(end)*Params.wY(1));
yactual = dot(Params.distpp(1:end-1),Params.ypp(1:end-1));
Params.ypp =  Params.ypp * ytarget / yactual;
Params.ypp(end) = ytop;

y = log(Params.ypp(:));
pi = Pi4(nmed,:);
Ey = pi*y;
Vary = pi*(y-Ey).^2;
Stdy = sqrt(Vary);

Ymat = Params.ypp(:)*Params.xY(:)';
Pmat = Params.distpp(:)*Params.wY(:)';
[Yqu,Yshares] = quints(Ymat,Pmat);


resid = [Vary - 0.015;
         Params.distpp(end)*Params.wY(end)-0.01;
         Yshares(5) - 0.6139];
         %sum(Yshares(7:8)) - (0.1637+0.1476)];
         % Yshares(8) - 0.1476];





function Pi = make_pi(pSwitch,pSwitch2,n)
Pi = eye(n);
for i=1:n
  if(i==1)
    Pi(i,i+1) = Pi(i,i+1) + pSwitch;
    Pi(i,i) = Pi(i,i) - pSwitch;
  elseif(i==n)
    Pi(i,i-1) = Pi(i,i-1) + pSwitch;
    Pi(i,i) = Pi(i,i) - pSwitch;
  elseif(i==n-1)
    Pi(i,i+1) = Pi(i,i+1) + pSwitch2;
    Pi(i,i-1) = Pi(i,i-1) + pSwitch;
    Pi(i,i) = Pi(i,i) - pSwitch - pSwitch2;
  else
    Pi(i,i+1) = Pi(i,i+1) + pSwitch;
    Pi(i,i-1) = Pi(i,i-1) + pSwitch;
    Pi(i,i) = Pi(i,i) - 2*pSwitch;
  end
end


function [resid,w] =  inciid_z(x)
global Params;
if(x(1)>=x(2))
  resid = 1e100;
end
if(Params.npp==1)
  p = [0.5 0.5];
else
  % p = [0.06 0.94];
  p = [0.25 0.75];
  p = [0.5 0.5];
end
logx = log(x);
Elogx = p*logx;
sigY = sqrt(0.061/Params.freq);
resid = [p*x - 1;0];
if(Params.npp==1)
  resid(2) = p*(logx-Elogx).^2 - sigY;
else
  resid(2) = x(1) - 0.2;
end
w = p';