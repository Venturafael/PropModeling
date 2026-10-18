function paths = propeller_repo_paths()
%PROPELLER_REPO_PATHS Localiza os dados sem depender da pasta atual.

utilities_dir = fileparts(mfilename('fullpath'));
matlab_root = fileparts(utilities_dir);

paths.matlab_root = matlab_root;
paths.functions = utilities_dir;
paths.data_root = fullfile(matlab_root, 'Data');
paths.apc_geometry = fullfile(paths.data_root, 'APC', 'geometry');
paths.apc_performance = fullfile(paths.data_root, 'APC', 'performance');
paths.geometry_22x10 = fullfile(paths.apc_geometry, '22x10-PERF.PE0');
paths.performance_22x10 = fullfile(paths.apc_performance, ...
    'PERFILES2', 'PER3_22x10.dat');
paths.polar_clark_y = fullfile(paths.data_root, 'ClarkY', ...
    'polares_clarky_expan.mat');
paths.airfoil_clark_y = fullfile(paths.data_root, 'ClarkY', 'clark_y.dat');
paths.matlab_results = fullfile(matlab_root, 'Results');
end
