function split_tissue_maps(varargin)

[compThreads, count] = sscanf(getenv('NSLOTS'), '%d');
if count == 1
  fprintf('NSLOTS=%d\n', compThreads);
  warning('off', 'MATLAB:maxNumCompThreads:Deprecated');
  autoCompThreads = maxNumCompThreads(compThreads);
end

if ~isdeployed
  addpath('./NIfTI_20140122');
  addpath(genpath('./TVAL3_v1.0'));
  addpath(genpath('./TVAL3D'));
end

SW_VER='1.0';

nargin = numel(varargin);
if nargin == 0
  fprintf('==========================================================================\n');
  fprintf('split_tissue_maps\n');
  fprintf('  Version %s\n', SW_VER);
  %fprintf('  Revision %s\n', SW_REV);
  fprintf('  Developed by Dongjin Kwon\n');
  fprintf('==========================================================================\n');
  fprintf('\n');
  fprintf('Required:\n');
  fprintf('  -i <file> : Input image.\n');
  fprintf('  -p <file> : Tissue probabilities. More than one need to be specified, each one preceded by its own -p.\n');
  fprintf('  -o <path> : Output prefix.\n');
  fprintf('Optional:\n');
  fprintf('  -mu <real> : Regualarization parameter (default: 10).\n');
  exit(0);
end

while numel(varargin) > 0 && isempty(varargin{1}),
  varargin(1) = [];
end
nargin = numel(varargin);

params.mu = 10;
p = 0;
a = 1;
while a <= nargin,
  switch varargin{a},
    case '-i',
      params.i = varargin{a+1};
      a = a + 2;
    case '-p',
      params.p{p+1} = varargin{a+1};
      p = p + 1;
      a = a + 2;
    case '-o',
      params.o = varargin{a+1};
      a = a + 2;
    case '-mu',
      params.mu = varargin{a+1};
      if ischar(params.mu),
        params.mu = str2double(params.mu);
      end
      a = a + 2;
    otherwise
      error('Unknown option: "%s"',varargin{a});
  end
end

fprintf('Input Image: %s\n', params.i);
for pi=1:p
  fprintf('Tissue Probability %d: %s\n', pi, params.p{pi});
end
fprintf('Output Prefix: %s\n', params.o);
fprintf('mu: %d\n', params.mu);

input_nii = load_untouch_nii(params.i);
for pi=1:p
  pbmap_nii{pi} = load_untouch_nii(params.p{pi});
end

dim_x = size(input_nii.img, 1);
dim_y = size(input_nii.img, 2);
dim_z = size(input_nii.img, 3);
n = dim_x * dim_y * dim_z;
nz = n * p;

A_i = zeros(nz, 1);
A_j = zeros(nz, 1);
A_v = zeros(nz, 1);
b = zeros(n, 1);

nz_idx = 1;
n_idx = 1;
for k=1:dim_z
  for j=1:dim_y
    for i=1:dim_x
      for pi=1:p
        A_i(nz_idx) = n_idx;
        A_j(nz_idx) = n * (pi-1) + n_idx;
        A_v(nz_idx) = pbmap_nii{pi}.img(i, j, k);
        nz_idx = nz_idx + 1;
      end

      b(n_idx) = input_nii.img(i, j, k);

      n_idx = n_idx + 1;
    end
  end
end

A_i(nz_idx:end) = [];
A_j(nz_idx:end) = [];
A_v(nz_idx:end) = [];

A = sparse(A_i, A_j, A_v, n, n * p, nz);

%alpha = A \ b;

opts.mu = params.mu;
opts.beta = 2^5;
opts.mu0 = opts.mu;
opts.beta0 = opts.beta;
opts.tol = 1e-6;
opts.tol_inn = 1e-3;
opts.maxit = 300;
opts.maxcnt = 10;
opts.TVnorm = 2;
opts.nonneg = true;
%opts.nonneg = false;
opts.TVL2 = true;
opts.isreal = false;
opts.scale_A = true;
opts.scale_b = true;
opts.disp = true;
opts.init = 1;

[alpha, out] = ftvcs_al_TVL2_3D(A, b, dim_x, dim_y, dim_z * p, opts);

for pi=1:p
  %out_file{pi} = strcat(params.o, '_', num2str(pi), '.nii.gz'); 
  out_file{pi} = strcat(params.o, '_mu', num2str(params.mu), '_', num2str(pi), '.nii.gz'); 
end

for pi=1:p
  alpha_nii = load_untouch_nii(params.i);
  alpha_nii.ext = [];

  alpha_nii.fileprefix = out_file{pi};
  for k=1:dim_z
    for j=1:dim_y
      for i=1:dim_x
        alpha_nii.img(i, j, k) = alpha(i, j, k + (pi-1) * dim_z);
      end
    end
  end
  
  save_untouch_nii(alpha_nii, out_file{pi});
end

exit(0);
