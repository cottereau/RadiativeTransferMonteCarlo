function config_h5_file = export_rtmc_config(geometry, source, material, observation, config_h5_file)
%% EXPORT_RTMC_CONFIG
% Exporta um objeto MaterialClass (já calculado ou não) e os parâmetros de simulação
% (geometry, source, observation) para um arquivo HDF5 (.h5) pronto para o código C++.
%
% Uso:
%   export_rtmc_config(geometry, source, material, observation);
%   export_rtmc_config(geometry, source, material, observation, 'meu_config.h5');
%
% Entradas:
%   - geometry: struct com 'dimension', 'frame', 'bnd'
%   - source  : struct com 'numberParticles', 'position', 'lambda'
%   - material: objeto MaterialClass preparado
%   - observation: struct com 'x', 'y', 'z', 'directions', 'time'
%   - config_h5_file (opcional): caminho para salvar o arquivo .h5 (default: 'cpp_config.h5')

if nargin < 5 || isempty(config_h5_file)
    config_h5_file = fullfile(pwd, 'cpp_config.h5');
end

fprintf('========================================================================\n');
fprintf(' EXPORTANDO CONFIGURAÇÃO MATERIAIS / SIMULAÇÃO PARA C++ (HDF5)\n');
fprintf(' Arquivo de saída: %s\n', config_h5_file);
fprintf('========================================================================\n');

% 1. Garantir que o MaterialClass está com todas as seções chocadas e preparadas
d = geometry.dimension;
if ~isprop(material, 'meanFreeTime') || isempty(material.meanFreeTime) || isempty(material.invcdf)
    fprintf('Preparando seções de choque (MaterialClass.prepareSigma)...\n');
    material = MaterialClass.prepareSigma(material, d);
end

% Remover arquivo antigo se existir
if exist(config_h5_file, 'file')
    delete(config_h5_file);
end

% 2. Geometria
h5create(config_h5_file, '/config/dimension', 1, 'Datatype', 'int32');
h5write(config_h5_file, '/config/dimension', int32(d));

frame_type = 2; % default esférico
if isfield(geometry, 'frame')
    if strcmpi(geometry.frame, 'cartesian'), frame_type = 0;
    elseif strcmpi(geometry.frame, 'cylindrical'), frame_type = 1;
    elseif strcmpi(geometry.frame, 'spherical'), frame_type = 2;
    end
end
h5create(config_h5_file, '/config/frame', 1, 'Datatype', 'int32');
h5write(config_h5_file, '/config/frame', int32(frame_type));

% Fronteiras (bnd)
bnd_dirs = []; bnd_vals = []; bnd_types = [];
if isfield(geometry, 'bnd') && ~isempty(geometry.bnd)
    for i = 1:length(geometry.bnd)
        bnd_dirs(end+1) = double(geometry.bnd(i).dir);
        bnd_vals(end+1) = double(geometry.bnd(i).val);
        type_val = 0; % 0 = refletiva
        if isfield(geometry.bnd(i), 'type') && strcmpi(geometry.bnd(i).type, 'absorbing')
            type_val = 1;
        end
        bnd_types(end+1) = type_val;
    end
end

if ~isempty(bnd_dirs)
    h5create(config_h5_file, '/config/bnd_dirs', length(bnd_dirs));
    h5write(config_h5_file, '/config/bnd_dirs', bnd_dirs(:));
    h5create(config_h5_file, '/config/bnd_vals', length(bnd_vals));
    h5write(config_h5_file, '/config/bnd_vals', bnd_vals(:));
    h5create(config_h5_file, '/config/bnd_types', length(bnd_types));
    h5write(config_h5_file, '/config/bnd_types', bnd_types(:));
end

% 3. Fonte
Nparticles = source.numberParticles;
h5create(config_h5_file, '/config/numberParticles', 1);
h5write(config_h5_file, '/config/numberParticles', double(Nparticles));

lambda_val = 2e-4;
if isfield(source, 'lambda') && ~isempty(source.lambda), lambda_val = source.lambda; end
h5create(config_h5_file, '/config/source_lambda', 1);
h5write(config_h5_file, '/config/source_lambda', double(lambda_val));

pos = [0 0 0];
if isfield(source, 'position')
    p = source.position;
    if length(p) == 2, pos = [p 0]; else pos = p; end
end
h5create(config_h5_file, '/config/source_position', 3);
h5write(config_h5_file, '/config/source_position', pos(:));

% Direction type: 0 = uniform, 1 = outgoing, 2 = upper
direction_type = 0;
if isfield(source, 'direction')
    if ischar(source.direction) || isstring(source.direction)
        if strcmpi(source.direction, 'outgoing'), direction_type = 1;
        elseif strcmpi(source.direction, 'upper'), direction_type = 2;
        end
    end
end
h5create(config_h5_file, '/config/source_direction_type', 1, 'Datatype', 'int32');
h5write(config_h5_file, '/config/source_direction_type', int32(direction_type));

% Polarization: 1 = P wave, 0 = S wave
polarization_val = 1;
if isfield(source, 'polarization')
    if strcmpi(source.polarization, 'S'), polarization_val = 0; end
end
h5create(config_h5_file, '/config/source_polarization', 1, 'Datatype', 'int32');
h5write(config_h5_file, '/config/source_polarization', int32(polarization_val));

% Radial source profile (invcdf)
if isfield(source, 'radial') && ~isempty(source.radial)
    Rmax = source.radial.GridVectors{1}(end);
    invcdfsource = compute_inverseCDF(source.radial, d, Rmax);
    x_cdf = invcdfsource.GridVectors{1}(:);
    y_r = invcdfsource.Values(:);

    h5create(config_h5_file, '/config/source_has_radial', 1, 'Datatype', 'int32');
    h5write(config_h5_file, '/config/source_has_radial', int32(1));

    h5create(config_h5_file, '/config/source_invcdf_x', length(x_cdf));
    h5write(config_h5_file, '/config/source_invcdf_x', double(x_cdf));

    h5create(config_h5_file, '/config/source_invcdf_y', length(y_r));
    h5write(config_h5_file, '/config/source_invcdf_y', double(y_r));
else
    h5create(config_h5_file, '/config/source_has_radial', 1, 'Datatype', 'int32');
    h5write(config_h5_file, '/config/source_has_radial', int32(0));
end

% 4. Material e Física
is_acoustics = material.acoustics;
h5create(config_h5_file, '/config/acoustics', 1, 'Datatype', 'int32');
h5write(config_h5_file, '/config/acoustics', int32(is_acoustics));

vp = 2.0; vs = 2.0;
if is_acoustics
    vp = material.v; vs = material.v;
else
    vp = material.vp; vs = material.vs;
end
h5create(config_h5_file, '/config/vp', 1); h5write(config_h5_file, '/config/vp', vp);
h5create(config_h5_file, '/config/vs', 1); h5write(config_h5_file, '/config/vs', vs);
h5create(config_h5_file, '/config/v', 1);  h5write(config_h5_file, '/config/v', vp);

% Fator de Qualidade Q & Frequência
QP = 1e12; QS = 1e12; freq = 10.0;
if isprop(material, 'Q') && ~isempty(material.Q)
    qval = material.Q(:);
    QP = qval(1);
    if length(qval) > 1, QS = qval(2); else QS = qval(1); end
end
if isprop(material, 'Frequency') && ~isempty(material.Frequency)
    freq = material.Frequency(1);
end

h5create(config_h5_file, '/config/QP', 1); h5write(config_h5_file, '/config/QP', double(QP));
h5create(config_h5_file, '/config/QS', 1); h5write(config_h5_file, '/config/QS', double(QS));
h5create(config_h5_file, '/config/frequency', 1); h5write(config_h5_file, '/config/frequency', double(freq));

timeSteps = 0;
if isprop(material, 'timeSteps') && ~isempty(material.timeSteps)
    timeSteps = material.timeSteps;
end
h5create(config_h5_file, '/config/timeSteps', 1, 'Datatype', 'int32');
h5write(config_h5_file, '/config/timeSteps', int32(timeSteps));

meanFreeTimeP = material.meanFreeTime(1);
meanFreeTimeS = meanFreeTimeP;
if length(material.meanFreeTime) > 1, meanFreeTimeS = material.meanFreeTime(2); end

h5create(config_h5_file, '/config/meanFreeTimeP', 1); h5write(config_h5_file, '/config/meanFreeTimeP', meanFreeTimeP);
h5create(config_h5_file, '/config/meanFreeTimeS', 1); h5write(config_h5_file, '/config/meanFreeTimeS', meanFreeTimeS);

P2P = 0.5; S2S = 0.5;
if isprop(material, 'P2P') && ~isempty(material.P2P), P2P = material.P2P; end
if isprop(material, 'S2S') && ~isempty(material.S2S), S2S = material.S2S; end
h5create(config_h5_file, '/config/P2P', 1); h5write(config_h5_file, '/config/P2P', P2P);
h5create(config_h5_file, '/config/S2S', 1); h5write(config_h5_file, '/config/S2S', S2S);

% Exportar tabelas CDF Inversa de Espalhamento (Funções de Fase)
if iscell(material.invcdf)
    export_invcdf_table(config_h5_file, '/config/invcdf_PP', material.invcdf{1,1});
    export_invcdf_table(config_h5_file, '/config/invcdf_PS', material.invcdf{1,2});
    export_invcdf_table(config_h5_file, '/config/invcdf_SP', material.invcdf{2,1});
    export_invcdf_table(config_h5_file, '/config/invcdf_SS', material.invcdf{2,2});
else
    export_invcdf_table(config_h5_file, '/config/invcdf_PP', material.invcdf);
end

% Exportar tabelas de reflexão de Zoeppritz
if ~is_acoustics
    try
        [out_Z, ~] = MaterialClass.ZoeppritzVsAir(material);
        x_deg = out_Z.j1_deg(:) / 90.0;
        y_rpp = out_Z.E_Rpp(:);
        y_rss = out_Z.E_Rsvsv(:);

        h5create(config_h5_file, '/config/Rpp_table_x', length(x_deg)); h5write(config_h5_file, '/config/Rpp_table_x', x_deg(:));
        h5create(config_h5_file, '/config/Rpp_table_y', length(y_rpp)); h5write(config_h5_file, '/config/Rpp_table_y', y_rpp(:));
        h5create(config_h5_file, '/config/Rss_table_x', length(x_deg)); h5write(config_h5_file, '/config/Rss_table_x', x_deg(:));
        h5create(config_h5_file, '/config/Rss_table_y', length(y_rss)); h5write(config_h5_file, '/config/Rss_table_y', y_rss(:));
    catch
        % Opcional
    end
end

% 5. Malha de Observação
x_bins = observation.x;
y_bins = [-pi pi]; if isfield(observation, 'y') && ~isempty(observation.y), y_bins = observation.y; end
z_bins = [-pi/2 pi/2]; if isfield(observation, 'z') && ~isempty(observation.z), z_bins = observation.z; end
dir_bins = [0 pi]; if isfield(observation, 'directions') && ~isempty(observation.directions), dir_bins = observation.directions; end
time_grid = observation.time;

h5create(config_h5_file, '/config/x_bins', length(x_bins)); h5write(config_h5_file, '/config/x_bins', x_bins(:));
h5create(config_h5_file, '/config/y_bins', length(y_bins)); h5write(config_h5_file, '/config/y_bins', y_bins(:));
h5create(config_h5_file, '/config/z_bins', length(z_bins)); h5write(config_h5_file, '/config/z_bins', z_bins(:));
h5create(config_h5_file, '/config/dir_bins', length(dir_bins)); h5write(config_h5_file, '/config/dir_bins', dir_bins(:));
h5create(config_h5_file, '/config/time_grid', length(time_grid)); h5write(config_h5_file, '/config/time_grid', time_grid(:));

fprintf('Exportação concluída com sucesso!\n');
fprintf('  - Dimensão : %dD\n', d);
fprintf('  - Tipo     : %s\n', ifelse(is_acoustics, 'Acústico', 'Elástico Bifásico'));
fprintf('  - Partículas: %e\n', Nparticles);
fprintf('  - Fronteiras: %d\n', length(bnd_dirs));
fprintf('Para rodar diretamente via terminal:\n');
fprintf('  ./cpp_rtmclib/build/rtmc_cpp "%s" "output.h5"\n', config_h5_file);
fprintf('========================================================================\n\n');

end

function export_invcdf_table(h5_file, base_path, invcdf_obj)
    if isempty(invcdf_obj) || ~isprop(invcdf_obj, 'GridVectors')
        x_grid = [0 1]; y_grid = [0 pi];
    else
        x_grid = invcdf_obj.GridVectors{1}(:);
        y_grid = invcdf_obj.Values(:);
    end
    h5create(h5_file, [base_path '_x'], length(x_grid));
    h5write(h5_file, [base_path '_x'], x_grid);
    h5create(h5_file, [base_path '_y'], length(y_grid));
    h5write(h5_file, [base_path '_y'], y_grid);
end

function out = ifelse(cond, a, b)
    if cond, out = a; else, out = b; end
end

function invcdf = compute_inverseCDF(f, d, Rmax)
    Nth = 10000;
    xth = linspace(0, Rmax, Nth);
    intF = trapz(xth, (xth.^(d-1)) .* f(xth));
    pdf = (xth.^(d-1)) .* f(xth) / intF;
    cdf = cumsum(pdf) * mean(diff(xth));
    cdf(1) = 0;
    cdf(end) = 1;
    ind = find(diff(cdf) > 0);
    ind = unique([ind ind+1]);

    [unique_cdf, idx] = unique(sort(cdf(ind)));
    ds = sort(xth(ind));
    unique_xth = ds(idx);

    if ~issorted(unique_cdf)
        [unique_cdf, sort_idx] = sort(unique_cdf);
        unique_xth = unique_xth(sort_idx);
    end

    invcdf = griddedInterpolant(unique_cdf, unique_xth, 'linear', 'nearest');
end
