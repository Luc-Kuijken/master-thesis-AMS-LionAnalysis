function mof = load_mof_structure(filepath, opts)
%LOAD_MOF_STRUCTURE Load MOF structure from .xyz file with unit cell vectors
%  mof = io.load_mof_structure(filepath, opts)
%  Returns struct with:
%    .name       : string with name of the MOF
%    .atoms      : cell array of atom types
%    .xyz        : Nby3 matrix of Cartesian coordinates [Angstrom]
%    .n_atoms    : number of atoms
%    .unit_cell  : 3by3 matrix of unit cell vectors [Angstrom]
%                  (VEC1, VEC2, VEC3 as rows)

    arguments
        filepath (1,1) string
        opts.Verbose (1,1) logical = true
    end
    
    if ~isfile(filepath)
        error("MOF file not found: %s", filepath); 
    end
    
    [~,mof_name,suffix] = fileparts(filepath);
    if opts.Verbose
        fprintf("Loading: %s%s\n", [mof_name,suffix])
    end
    
    fid = fopen(filepath, 'r');
    
    % Read number of atoms
    n_atoms = str2double(fgetl(fid));
    
    % Skip comment line
    fgetl(fid);
    
    % Preallocate
    atoms = cell(n_atoms, 1);
    xyz = zeros(n_atoms, 3);
    
    % Read atomic data
    for i = 1:n_atoms
        line = fgetl(fid);
        parts = strsplit(strtrim(line));
        atoms{i} = parts{1};
        xyz(i,:) = [str2double(parts{2}), str2double(parts{3}), str2double(parts{4})];
    end
    
    % Initialize unit cell (default to max values of coords)
    maximal_coords = max(xyz);
    unit_cell = eye(3) .* maximal_coords;
    
    % Read remaining lines for unit cell vectors
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line), break; end
        line = strtrim(line);
        
        if startsWith(line, 'VEC1')
            parts = strsplit(line);
            unit_cell(1,:) = [str2double(parts{2}), str2double(parts{3}), str2double(parts{4})];
        elseif startsWith(line, 'VEC2')
            parts = strsplit(line);
            unit_cell(2,:) = [str2double(parts{2}), str2double(parts{3}), str2double(parts{4})];
        elseif startsWith(line, 'VEC3')
            parts = strsplit(line);
            unit_cell(3,:) = [str2double(parts{2}), str2double(parts{3}), str2double(parts{4})];
        end
    end
    
    fclose(fid);
    
    % Store in struct
    mof.name = mof_name;
    mof.atoms = atoms;
    mof.xyz = xyz;
    mof.n_atoms = n_atoms;
    mof.unit_cell = unit_cell;
end