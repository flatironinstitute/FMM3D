%testing FMM library
clear srcinfo1 

disp("------ Testing hfmm -------");
FileName1 = '../data/CombinedMesh_bary_101309.mat';
% FileName1 = 'CombinedMesh_bary2_101309.mat';
[P1, t1, normals1, Area1, Center1, Indicator1, Size1, remove1] = MyParLoad(FileName1);

CC = Center1';

nmax = 8500000;
nsources = min(length(CC(1,:)),nmax);

CC = CC(:,1:nsources);

xmax = max(CC(1,:));
xmin = min(CC(1,:));


ymax = max(CC(2,:));
ymin = min(CC(2,:));


zmax = max(CC(3,:));
zmin = min(CC(3,:));

bx = xmax-xmin; by = ymax-ymin; bz = zmax-zmin;

bsize = max(bx,by);
bsize = max(bz,bsize);

kb = 100;

srcinfo1.sources = CC;
charges		 = rand(1,nsources)*10^-18;
srcinfo1.charges = charges;

k = kb/bsize+0.0000000001j; % wavenumber

% 2 digits accuracy new
% potentai plus field
eps = 0.5e-2;
pg = 2;	 
tic
H1 = hfmm3d(eps,k,srcinfo1,pg);
newFMM_time = toc