function plot_mof_frame(mof, opts)
%PLOT_MOF_FRAME Plot MOF structure frame (unit cell box) on current axes
%  plotting.plot_mof_frame(mof, Name=Value)
%  Inputs: mof struct from io.load_mof_structure
%  Name-Value options:
%    Projection  : "XY" | "YZ" | "ZX" | "3D" (default: "3D")
%    Color       : line color for unit cell box (default: [0.3 0.3 0.3] dark gray)
%    LineWidth   : line width for unit cell box (default: 1.5)
%    LineStyle   : line style for unit cell box (default: '-')
%    AtomScale   : scaling factor for atom sizes (default: 1)
%    AtomAlpha   : opacity value for atoms (default: 1)

    arguments
        mof (1,1) struct
        opts.Projection (1,1) string = "3D"
        opts.Color (1,3) double = [0.3 0.3 0.3]
        opts.LineWidth (1,1) double = 1
        opts.LineStyle (1,1) string = "--"
        opts.AtomScale (1,1) double = 1
        opts.AtomAlpha (1,1) double = 1
    end
    
    % Get unit cell vectors
    a = mof.unit_cell(1,:);  % VEC1
    b = mof.unit_cell(2,:);  % VEC2
    c = mof.unit_cell(3,:);  % VEC3
    
    % Define unit cell vertices (8 corners)
    origin = [0, 0, 0];
    vertices = [
        origin;           % 1: (0,0,0)
        a;                % 2: (a,0,0)
        a + b;            % 3: (a,b,0)
        b;                % 4: (0,b,0)
        c;                % 5: (0,0,c)
        a + c;            % 6: (a,0,c)
        a + b + c;        % 7: (a,b,c)
        b + c             % 8: (0,b,c)
    ];
    
    % Define ALL edges (for 3D or reference)
    all_edges = [
        1,2; 2,3; 3,4; 4,1;  % bottom face
        5,6; 6,7; 7,8; 8,5;  % top face
        1,5; 2,6; 3,7; 4,8   % vertical edges
    ];
    
    % Select edges and dimensions based on projection
    switch opts.Projection
        case "XY"
            % Project onto XY plane: show outer rectangle
            edges = [1,2; 2,3; 3,4; 4,1];  % bottom face only
            dim1 = 1; dim2 = 2;  % X, Y
        case "YZ"
            % Project onto YZ plane: show outer rectangle
            edges = [4,1; 1,5; 5,8; 8,4];  % left face
            dim1 = 2; dim2 = 3;  % Y, Z
        case "ZX"
            % Project onto ZX plane: show outer rectangle
            edges = [1,2; 2,6; 6,5; 5,1];  % front face
            dim1 = 3; dim2 = 1;  % Z, X
        case "3D"
            % Full 3D: use all edges
            edges = all_edges;
            dim1 = []; dim2 = [];
        otherwise
            error('Projection must be "XY", "YZ", "ZX", or "3D"');
    end
    
    % Plot box using unified helper
    if opts.Projection == "3D"
        plot_unitcell(vertices, edges, opts, true);
        
        % Optionally plot atoms in 3D
        if opts.AtomAlpha > 0
            plot_atoms_3d(mof, opts);
        end
    else
        plot_unitcell(vertices(:,[dim1, dim2]), edges, opts, false);
        
        % Optionally plot atoms in 2D
        if opts.AtomAlpha > 0
            plot_atoms_2d(mof, opts, dim1, dim2);
        end
    end
    
end

%% Helper: Unified box plotting (2D and 3D)
function plot_unitcell(vertices, edges, opts, is3d)
    %  vertices: Nx2 (2D) or Nx3 (3D) matrix
    %  edges: Mx2 matrix of vertex index pairs
    %  opts: struct with Color, LineStyle, LineWidth
    %  is3d: logical, true for 3D plotting

    for i = 1:size(edges, 1)
        idx1 = edges(i,1);
        idx2 = edges(i,2);
        v1 = vertices(idx1, :);
        v2 = vertices(idx2, :);
        
        if is3d
            plot3([v1(1), v2(1)], [v1(2), v2(2)], [v1(3), v2(3)], ...
                  'LineStyle', opts.LineStyle, 'Color', opts.Color, ...
                  'LineWidth', opts.LineWidth, 'HandleVisibility', 'off');
        else
            plot([v1(1), v2(1)], [v1(2), v2(2)], ...
                 'LineStyle', opts.LineStyle, 'Color', opts.Color, ...
                 'LineWidth', opts.LineWidth, 'HandleVisibility', 'off');
        end
    end
end

%% Helper: Get atom properties (color, size)
function [color, size, edge_color] = get_atom_properties(atom_type, base_scale)
%GET_ATOM_PROPERTIES Return color and size for specific atom type
%  [color, size] = get_atom_properties(atom_type, base_scale)
%  Atom-specific styling based on UiO-66 structure
%  Zr: yellow/gold, large
%  C:  gray, medium
%  O:  red, small
%  H:  not plotted (returns empty)

    atom_type = upper(strtrim(atom_type));
    
    switch atom_type
        case 'ZR'
            color = [0.9, 0.7, 0.1];    % Gold
            size = 10 * base_scale;     % Large
        case 'C'
            color = [0.5, 0.5, 0.5];    % Gray
            size = 5 * base_scale;     	% Medium
        case 'O'
            color = [0.9, 0.1, 0.1];    % Red
            size = 4 * base_scale;    	% Small
        case 'H'
            color = [];                	% Don't plot hydrogen
            size = 0;
        otherwise
            % Unknown atom: light purple, small
            color = [0.7, 0.5, 0.7];
            size = 3 * base_scale;
    end
    edge_color = color * 0.8;
end

%% Helper: Plot 2D atoms
function plot_atoms_2d(mof, opts, dim1, dim2)
    coords = mof.xyz;
    atoms = mof.atoms;
    
    % Plot each atom type separately for proper coloring
    for i = 1:numel(atoms)
        [color, size, edge_color] = get_atom_properties(atoms{i}, opts.AtomScale);
        
        % Skip hydrogen (or if size is 0)
        if isempty(color) || size == 0
            continue
        end
        
        scatter(coords(i, dim1), coords(i, dim2), size^2, ...
             'MarkerFaceColor', color, 'MarkerEdgeColor', edge_color,...
             'MarkerFaceAlpha', opts.AtomAlpha, 'MarkerEdgeAlpha', opts.AtomAlpha, ...
             'LineWidth', 0.5, 'HandleVisibility', 'off');
    end
end

%% Helper: Plot 3D atoms
function plot_atoms_3d(mof, opts)
    coords = mof.xyz;
    atoms = mof.atoms;
    
    % Plot each atom type separately for proper coloring
    for i = 1:numel(atoms)
        [color, size, edge_color] = get_atom_properties(atoms{i}, opts.AtomScale);
        
        % Skip hydrogen (or if size is 0)
        if isempty(color) || size == 0
            continue
        end
        
        scatter3(coords(i, 1), coords(i, 2), coords(i, 3), size^2, ...
             'MarkerFaceColor', color, 'MarkerEdgeColor', edge_color,...
             'MarkerFaceAlpha', opts.AtomAlpha, 'MarkerEdgeAlpha', opts.AtomAlpha, ...
             'LineWidth', 0.5, 'HandleVisibility', 'off');
    end
end