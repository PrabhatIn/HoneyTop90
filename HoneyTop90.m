function HoneyTop90(HNex,HNey,volfrac,penal,rfill,ft)
%% ---------- Material properties ----------
E0 = 1;
Emin = E0*1e-9;
%% ---Element connectivity, nodal coordinates, Finite element analysis preparation--
NstartVs = reshape(1:(1+2*HNex)*(1+HNey),1+2*HNex,1+HNey);
DOFstartVs = reshape(2*NstartVs(1:end-1,1:end-1)-1,2*HNex*HNey,1);
NodeDOFs = repmat(DOFstartVs,1,8) + repmat([2*(2*HNex+1) + [2 3 0 1] 0 1 2 3 ],2*HNex*HNey,1);
ActualDOFs = NodeDOFs(setdiff(1:2*HNex*HNey,(2*HNex:2*HNex:2*HNex*HNey)' +  mod(1:HNey,2)'),:);
HoneyDOFs = [ActualDOFs(2:2:end,1:2), ActualDOFs(1:2:end,:), ActualDOFs(2:2:end,7:8)];
Ncyi = repmat(reshape(repmat([-0.25 0.25]',HNey+1,1),2,HNey+1)+reshape(1.5*sort(repmat((0:HNey)',2,1)),2,(HNey+1)),HNex+1,1);
Ncyi(:,1:2:end) = flip(Ncyi(:,1:2:end));
Ncyf = Ncyi(1:end-1,:); % final arrays containing y-coordinates
HoneyNCO=(1/sqrt(3))*[repmat((0:cos(pi/6):2*HNex*cos(pi/6)),1,HNey+1)' Ncyf(:)];%node co
if(mod(HNey,2)==0)
 HoneyDOFs(end:-1:end-(HNex)+2,1:6) = HoneyDOFs(end:-1:end-(HNex)+2,1:6)-2; % Updating
 HoneyNCO([(2*HNex+1)*HNey+1;(2*HNex+1)*(HNey+1)],:) = []; % Removing hangining nodes
end
[Nelem,Nnode] = deal(size(HoneyDOFs,1),size(HoneyNCO,1)); % elem #, node #
 F = sparse(2*((2*HNex+1)*HNey+1),1,-1,2*Nnode,1);        % Applied load
 U = zeros(2*Nnode,1);                                    % Initializing U
 fixeddofs = [2*(1:2*HNex+1:(2*HNex+1)*HNey+1)-1,(2*(2*HNex+1))];  % Fixed DOFs
 alldofs = 1:2*Nnode;                                              % Total DOFs
 freedofs = setdiff(alldofs,fixeddofs);                            % Free DOFs
 iK = reshape(kron(HoneyDOFs,ones(12,1))',144*Nelem,1);
 jK = reshape(kron(HoneyDOFs,ones(1,12))',144*Nelem,1);
 KE = E0*[616.43012 92.77147 -168.07333 65.54377 -232.28511 -0.00032 -120.65312 -83.28564 -71.60020 -92.77115 -23.81836 17.74187;
 92.77147 509.30685 101.02751 -71.90335 0.00032 -18.03857 -83.28564 -24.48314 -92.77179 -178.72347 -17.74187 -216.15832;
-168.07333 101.02751 455.74522 0.00000 -168.07333 -101.02751 -71.60020 -92.77179 23.60185 -0.00000 -71.60020 92.77179;
65.54377 -71.90335 0.00000 669.99176 -65.54377 -71.90335 -92.77115 -178.72347 -0.00000 -168.73811 92.77115 -178.72347;
-232.28511 0.00032 -168.07333 -65.54377 616.43012 -92.77147 -23.81836 -17.74187 -71.60020 92.77115 -120.65312 83.28564;
-0.00032 -18.03857 -101.02751 -71.90335 -92.77147 509.30685 17.74187 -216.15832 92.77179 -178.72347 83.28564 -24.48314;
-120.65312 -83.28564 -71.60020 -92.77115 -23.81836 17.74187 616.43012 92.77147 -168.07333 65.54377 -232.28511 -0.00032;
-83.28564 -24.48314 -92.77179 -178.72347 -17.74187 -216.15832 92.77147 509.30685 101.02751 -71.90335 0.00032 -18.03857;
-71.60020 -92.77179 23.60185 -0.00000 -71.60020 92.77179 -168.07333 101.02751 455.74522 0.00000 -168.07333 -101.02751;
-92.77115 -178.72347 -0.00000 -168.73811 92.77115 -178.72347 65.54377 -71.90335 0.00000 669.99176 -65.54377 -71.90335;
-23.81836 -17.74187 -71.60020 92.77115 -120.65312 83.28564 -232.28511 0.00032 -168.07333 -65.54377 616.43012 -92.77147;
17.74187 -216.15832 92.77179 -178.72347 83.28564 -24.48314 -0.00032 -18.03857 -101.02751 -71.90335 -92.77147 509.30685;]/1000;% elem stiffness 
%% ---------- Filter preperation ---------
Cxx = repmat([sqrt(3)/2*(1:2:2*HNex-1) sqrt(3)*(1:1:HNex-1)],1,ceil(HNey/2));
Cyy= (repmat(3/4,HNex,HNey) + repmat(3/2*(0:HNey-1),HNex,1));
Cyy(HNex+1:2*HNex:length(Cyy(:))) = [];
ct = [Cxx(1:length(Cyy))' Cyy']*(1/sqrt(3));             % Centre coordinates
DD = cell(Nelem,1);                                      % Initializing 
for j = 1:Nelem
    Cent_dist = sqrt((ct(j,1)-ct(:,1)).^2+((ct(j,2)-ct(:,2)).^2));
    [I,J] = find(Cent_dist<=rfill);
    DD{j} = [I,J+(j-1),Cent_dist(I)];
end
 DD = cell2mat(DD);
 HHs = sparse(DD(:,2),DD(:,1),1-DD(:,3)/rfill);
 HHs = spdiags(1./sum(HHs,2),0,Nelem,Nelem)*HHs;                           
%% ---------- Initialization -------------
x = volfrac*ones(Nelem,1);                                               % Initial guess
[xPhys,loop,change,maxiter,dv,move] = deal(x,0,1,200,ones(Nelem,1),0.2); % Parameters 
%% ---------- Start optimization ----------
while (change > 0.01 && loop < maxiter)
  loop = loop + 1;
  %% ---------- Finite element analysis ----------
  sK = reshape(KE(:)*(Emin + xPhys'.^penal*(E0-Emin)),144*Nelem,1);
  K = sparse(iK,jK,sK);                                            % Global stiffness 
  U(freedofs) = decomposition(K(freedofs,freedofs),'chol','lower')\F(freedofs);
  %% ---------- Objective and sensitivities evaluation ----------
  ce = sum((U(HoneyDOFs)*KE).*U(HoneyDOFs),2);
  c =  sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));               % Finding objective 
  dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;                     % Obj. sensitivities
   %% ---------- Using Filters ---------------------------
   if ft == 1
    dc = HHs'*(x.*dc)./max(1e-3,x);
  elseif ft == 2
    dc = HHs'*(dc);
    dv = HHs'*(dv);
   end
  %% ---------- Optimality criteria update ----------------
   xOpt = x;
  [xUpp, xLow] = deal (xOpt + move, xOpt - move);              % Upp. & low. limits
  OcC = xOpt.*(sqrt(-dc./dv));                                 % Opt. parameter 
  inL = [0, mean(OcC)/volfrac];                                % Lag. Mul. range
  while (inL(2)-inL(1))/(inL(2)+inL(1))> 1e-3
    lmid = 0.5*(inL(2)+ inL(1));
    x = max(0,max(xLow,min(1,min(xUpp,OcC/lmid))));
    if mean(x)>volfrac, inL(1) = lmid; else, inL(2) = lmid; end 
  end
  if ft == 1 || ft ==0, xPhys = x; elseif ft == 2, xPhys = HHs'*x; end
  change = max(abs(xOpt-x));
  %% ---------- Results printing ----------------------------
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, mean(xPhys),change);
  %% ---------- Plotting intermediate designs ---------------
  colormap('gray'); scatter(ct(:,1),ct(:,2),[],1-xPhys,'filled'); axis off equal; pause(1e-6);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This Matlab code is presented for the educational purposes and its      %
%  details can be found in "HoneyTop90: A 90-line MATLAB code for          %
%  topology optimization using honeycomb tessellation" Optimization and    %
%  Engineering Journal, 2022, in press                                     %         
%  Please write : prabhatk@iisc.ac.in or prabhatkumar.rns@gmail.com        %
%  for any comments                                                        %
%                                                                          %                                              
%  Disclaimer:                                                             %
%  The author reserves all rights but does not guaranty that the code is   %
%  free from errors and the author shall not be liable in any event        %
%  caused by the use of the program.                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
