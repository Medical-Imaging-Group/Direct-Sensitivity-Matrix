
function [fwd_mesh,pj_error] = reconstruct_stnd_cw_region_direct_sensitivity(fwd_mesh,...
						  data_fn,...
						  iteration,...
						  lambda,...
						  output_fn,...
						  filter_n,...
						  region)

% [fwd_mesh,pj_error] = reconstruct_stnd_cw_region(fwd_mesh,...
%						  data_fn,...
%						  iteration,...
%						  lambda,...
%						  output_fn,...
%						  filter_n,...
%						  region)
%
% CW Reconstruction program for standard meshes using region based priori
% Needs meshes that have region labels
% Assumes that Intensity only data, and only reconstructs mua
%
% fwd_mesh is the input mesh (variable or filename)
% data_fn is the boundary data (variable or filename)
% iteration is the max number of iterations
% lambda is the initial regularization value
%    e.g. [20 20 15; 20 20 15] for a 3D mesh with 2 regions
% output_fn is the root output filename
% filter_n is the number of mean filters
% region is an array of the regions (e.g. [0 1 2])



% set modulation frequency to zero.
frequency = 0;

tic;

%****************************************
% If not a workspace variable, load mesh
if ischar(fwd_mesh)== 1
  fwd_mesh = load_mesh(fwd_mesh);
end
if ~strcmp(fwd_mesh.type,'stnd')
    errordlg('Mesh type is incorrect','NIRFAST Error');
    error('Mesh type is incorrect');
end

% read data - This is the calibrated experimental data or simulated data
anom = load_data(data_fn);
if ~isfield(anom,'paa')
    errordlg('Data not found or not properly formatted','NIRFAST Error');
    error('Data not found or not properly formatted');
end
anom = anom.paa;
% anom = log(anom(:,1));
anom = log(anom(:,1));
% find NaN in data

datanum = 0;
[ns,junk]=size(fwd_mesh.source.coord);
for i = 1 : ns
  for j = 1 : length(fwd_mesh.link(i,:))
      datanum = datanum + 1;
      if fwd_mesh.link(i,j) == 0
          anom(datanum,:) = NaN;
      end
  end
end

ind = find(isnan(anom(:,1))==1);
% set mesh linkfile not to calculate NaN pairs:
link = fwd_mesh.link';
link(ind) = 0;
fwd_mesh.link = link';
recon_mesh.link = link';
clear link
% remove NaN from data
ind = setdiff(1:size(anom,1),ind);
anom = anom(ind,:);
clear ind;

% Initiate projection error
pj_error = [];

% Initiate log file
fid_log = fopen([output_fn '.log'],'w');
fprintf(fid_log,'Forward Mesh   = %s\n',fwd_mesh.name);
if ischar(data_fn) ~= 0
    fprintf(fid_log,'Data File      = %s\n',data_fn);
end
fprintf(fid_log,'Initial Reg    = %d\n',lambda);
fprintf(fid_log,'Filter         = %d\n',filter_n);
fprintf(fid_log,'Output Files   = %s_mua.sol\n',output_fn);
fprintf(fid_log,'               = %s_mus.sol\n',output_fn);

% This calculates the mapping matrix that reduces Jacobian from nodal
% values to regional values
disp('calculating regions');
if ~exist('region','var')
    region = unique(fwd_mesh.region);
end
K = region_mapper(fwd_mesh,region);


for it = 1 : iteration
  
 
  data=femdata(fwd_mesh,0);
  clear ref;

 ref(:,1) = log(data.amplitude);

  data_diff = (anom-ref);

  pj_error = [pj_error sum(abs(data_diff.^2))]; 
  
  disp('---------------------------------');
  disp(['Iteration Number          = ' num2str(it)]);
  disp(['Projection error          = ' num2str(pj_error(end))]);

  fprintf(fid_log,'---------------------------------\n');
  fprintf(fid_log,'Iteration Number          = %d\n',it);
  fprintf(fid_log,'Projection error          = %f\n',pj_error(end));

  if it ~= 1
    p = (pj_error(end-1)-pj_error(end))*100/pj_error(end-1);
    disp(['Projection error change   = ' num2str(p) '%']);
    fprintf(fid_log,'Projection error change   = %f %%\n',p);
    if (p) <= 2
      disp('---------------------------------');
      disp('STOPPING CRITERIA REACHED');
      fprintf(fid_log,'---------------------------------\n');
      fprintf(fid_log,'STOPPING CRITERIA REACHED\n');
     break
    end
  end
   
  if it==1
   regions = unique(fwd_mesh.region);
NR = length(regions);
R = zeros(length(fwd_mesh.mua), NR);
MuaREC = zeros(NR,1);
MusREC = MuaREC;

% store the region indices
for i = 1:NR
    R1 = find(fwd_mesh.region == regions(i));
    RL(i) = length(R1);
    R(1:length(R1), i) = R1;
    MuaREC(i) = mean(fwd_mesh.mua(R1));
    MusREC(i) = mean(fwd_mesh.mus(R1));
    clear R1;
end
pert1 = 0.1/100;%0.1% Perturbation
dmua=zeros(1,length(NR));
t1=sum(K);
  end 
  




 N = fwd_mesh.mua'*K./t1;
 
%%% Estimation using Direct-Sensitivity Method-Matrix free estimation of mu_a%%%%%%%%%%%%%

for i=1:length(N)
ind = R(1:RL(i), i);    
p_mesh=fwd_mesh;

    delmua1 = (mean(p_mesh.mua(ind))*pert1);
    p_mesh.mua(ind)=p_mesh.mua(ind)+ delmua1;
    
    
    I_p=femdata(p_mesh,0);
    ln_p1 = log(I_p.amplitude);
    
    delp = (ln_p1-ref);
        Ind = abs(delp) < max(abs(delp));

%      Ind1 = find(abs(delp) >= max(abs(delp)));
%     Ind1
     delp(Ind) = 0;
    
    Jinv =((delmua1))./delp;
    Jinv(isinf(Jinv)) =0;
    
   
    N2 = 1./ ((Jinv./N(i))'*ref);

    
    dmua(i) = (Jinv'*data_diff).*N2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%



  % use region mapper to unregionize!
  foo = K*dmua';

  % Update values
  fwd_mesh.mua = fwd_mesh.mua + foo;
  fwd_mesh.kappa = 1./(3.*(fwd_mesh.mus + fwd_mesh.mua));
  
%%%%%%%%%%DSM Ends%%%%%%%%%%%%%
  
  clear foo 
  % We dont like -ve mua or mus! so if this happens, terminate
  if (any(fwd_mesh.mua<0) | any(fwd_mesh.mus<0))
    disp('---------------------------------');
    disp('-ve mua or mus calculated...not saving solution');
    fprintf(fid_log,'---------------------------------\n');
    fprintf(fid_log,'STOPPING CRITERIA REACHED\n');
    break
  end
  
  % Filtering if needed!
  if filter_n > 1
    fwd_mesh = mean_filter(fwd_mesh,abs(filter_n));
  elseif filter_n < 0
    fwd_mesh = median_filter(fwd_mesh,abs(filter_n));
  end

  if it == 1
    fid = fopen([output_fn '_mua.sol'],'w');
  else
    fid = fopen([output_fn '_mua.sol'],'a');
  end
  fprintf(fid,'solution %g ',it);
  fprintf(fid,'-size=%g ',length(fwd_mesh.nodes));
  fprintf(fid,'-components=1 ');
  fprintf(fid,'-type=nodal\n');
  fprintf(fid,'%f ',fwd_mesh.mua);
  fprintf(fid,'\n');
  fclose(fid);
  
  if it == 1
    fid = fopen([output_fn '_mus.sol'],'w');
  else
    fid = fopen([output_fn '_mus.sol'],'a');
  end
  fprintf(fid,'solution %g ',it);
  fprintf(fid,'-size=%g ',length(fwd_mesh.nodes));
  fprintf(fid,'-components=1 ');
  fprintf(fid,'-type=nodal\n');
  fprintf(fid,'%f ',fwd_mesh.mus);
  fprintf(fid,'\n');
  fclose(fid);
end

% close log file!
time = toc;
fprintf(fid_log,'Computation Time = %f\n',time);
disp(['Computation Time = ' num2str(time) 'sec']);
fclose(fid_log);




function KKK=region_mapper(mesh,region)

% created a matrix which is in effect a mapper to move 
% between nodal basis or region basis.
% h dehghani May 8 2002

%nregions = max(mesh.region)+1;
nregion = length(region);
disp(['Number of regions = ' num2str(nregion)]);
nnodes = length(mesh.nodes);

% create empty mapping matrix
K = sparse(nnodes,nregion);

% Assign mapping functions, for each node belonging to each region
for j = 1 : nregion
  K(find(mesh.region==region(j)),j) = 1;
end

% find the total number of assigned nodes
N = full(sum(sum(K)));

% Here if some node is not in region, must account for it
if N ~= length(mesh.nodes)
  KK = sparse(nnodes,nnodes-N);
  for k = 1 : length(region)
    if k == 1
      a = find(mesh.region~=region(k));
    else
      a = intersect(find(mesh.region~=region(k)),a);
    end
  end
  for i = 1 : length(a)
    KK(a(i),i) = 1;
  end
  KKK = [K KK];
else
  KKK = K;
end
