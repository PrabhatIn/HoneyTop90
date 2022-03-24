function HoneyGen(HNex,HNey,a)
%% Example: generating mesh grid with a (edge-lenght) = 1 unit
%HoneyGen(100,100,1)
%% ---Element connectivity, nodal coordinates--
NstartVs = reshape(1:(1+2*HNex)*(1+HNey),1+2*HNex,1+HNey);
DOFstartVs = reshape(2*NstartVs(1:end-1,1:end-1)-1,2*HNex*HNey,1);
NodeDOFs = repmat(DOFstartVs,1,8) + repmat([2*(2*HNex+1) + [2 3 0 1] 0 1 2 3 ],2*HNex*HNey,1);
ActualDOFs = NodeDOFs(setdiff(1:2*HNex*HNey,(2*HNex:2*HNex:2*HNex*HNey)' +  mod(1:HNey,2)'),:);
HoneyDOFs = [ActualDOFs(2:2:end,1:2), ActualDOFs(1:2:end,:), ActualDOFs(2:2:end,7:8)];
Ncyi = repmat(reshape(repmat([-0.25 0.25]',HNey+1,1),2,HNey+1)+reshape(1.5*sort(repmat((0:HNey)',2,1)),2,(HNey+1)),HNex+1,1);
Ncyi(:,1:2:end) = flip(Ncyi(:,1:2:end));
Ncyf = Ncyi(1:end-1,:); % final arrays containing y-coordinates
HoneyNCO = a*[repmat((0:cos(pi/6):2*HNex*cos(pi/6)),1,HNey+1)' Ncyf(:)];%node co
if(mod(HNey,2)==0)
 HoneyDOFs(end:-1:end-(HNex)+2,1:6) = HoneyDOFs(end:-1:end-(HNex)+2,1:6)-2; % Updating
 HoneyNCO([(2*HNex+1)*HNey+1;(2*HNex+1)*(HNey+1)],:) = []; % Removing hangining nodes
end
HoneyElem = HoneyDOFs(:,2:2:end)/2;% element connectivity matrix
%% plotting the mesh
  patch('Faces',HoneyElem,'Vertices',HoneyNCO,'FaceColor','white','EdgeColor','k'); axis off equal; pause(1e-6);
 %% Finding the total number of elements and nodes in the mesh grid
 %[Nelem,Nnode] = deal(size(HoneyDOFs,1),size(HoneyNCO,1)); % elem #, node #
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This Matlab code is presented for the educational purposes and its      %
%  details can be found in "HoneyTop90: A 90-line MATLAB code for          %
%  topology optimization using honeycomb tessellation"  Optimization and   %
%  Engineering Journal, 2022, in press                                     %
%  Please write : prabhatk@iisc.ac.in or prabhatkumar.rns@gmail.com        %
%  for any comments                                                        %
%                                                                          %                                              
%  Disclaimer:                                                             %
%  The author reserves all rights but does not guaranty that the code is   %
%  free from errors and the author shall not be liable in any event        %
%  caused by the use of the program.                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
