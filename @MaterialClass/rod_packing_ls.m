function [centers, num_rods] = rod_packing_ls(L, DD, target_phi)
% rod_packing_ls_1D.m
% 1D version of RSA + Lubachevsky-Stillinger grow-and-relax algorithm
% for hard rod packing on a periodic line.

t_start = tic;

% =========================================================
% 1. Parameters
% =========================================================
escala = 1.0;
% L        = 1.0 * escala;      % Box length
% DD       = 0.02 * escala;     % Target rod length (diameter)
% target_phi = 0.45;

vol_rod   = DD;                     % "volume" (length) of one rod
vol_box   = L;
num_rods  = ceil((vol_box * target_phi) / vol_rod);

fprintf('=============================================\n');
fprintf(' 1D Rod Packing: RSA + LS Compression\n');
fprintf('=============================================\n');
fprintf(' Target phi      : %g\n', target_phi);
fprintf(' Target N rods   : %d\n', num_rods);
fprintf(' Box length      : %g\n', L);
fprintf(' Final length    : %g\n', DD);

centers = zeros(num_rods, 1);   % Column vector for 1D positions

rng('shuffle');

% =========================================================
% 2. Phase 1 — RSA at reduced diameter
% =========================================================
DD_current = DD * (0.42 / target_phi);   % 1D scaling: linear (exponent = 1)

fprintf('\n--- Phase 1: RSA at reduced length ---\n');
fprintf('    RSA length     : %g\n', DD_current);
fprintf('    RSA phi target : %g\n', (num_rods * DD_current) / L);

[cell_head, cell_next, cell_size, ncell] = build_cell_list_1D(DD_current, L, num_rods);

max_att  = 500000000;
attempts = 0;
cnt      = 0;

while cnt < num_rods && attempts < max_att
    cand = rand * L;
    
    overlapping = false;
    ic = min(floor(cand / cell_size) + 1, ncell);
    
    for dx = -1:1
        jc = mod(ic + dx - 1, ncell) + 1;
        ni = cell_head(jc);
        
        while ni ~= 0
            d_vec = abs(cand - centers(ni));
            if d_vec > L*0.5
                d_vec = L - d_vec;
            end
            if d_vec < DD_current
                overlapping = true;
                break;
            end
            ni = cell_next(ni);
        end
        if overlapping, break; end
    end

    if ~overlapping
        cnt = cnt + 1;
        centers(cnt) = cand;
        
        ic = min(floor(cand / cell_size) + 1, ncell);
        cell_next(cnt) = cell_head(ic);
        cell_head(ic) = cnt;
        
        if mod(cnt, 1000) == 0
            fprintf('  RSA placed: %d / %d\n', cnt, num_rods);
        end
        attempts = 0;
    else
        attempts = attempts + 1;
    end
end

if cnt < num_rods
    error('ERROR: RSA could not place all rods at reduced length!');
end
fprintf('  RSA complete: %d rods placed.\n', cnt);

% =========================================================
% 3. Phase 2 — Grow lengths + MC displacement to relax
% =========================================================
fprintf('\n--- Phase 2: Grow & Relax (LS-MC) ---\n');

grow_rate   = 1.005;
max_iter_mc = 5000;
step_size   = DD * 0.05;

while DD_current < DD
    DD_current = min(DD_current * grow_rate, DD);

    [cell_head, cell_next, cell_size, ncell] = build_cell_list_1D(DD_current, L, num_rods);
    [cell_head, cell_next] = rebuild_cell_list_from_centers_1D(centers, num_rods, cell_size, ncell);

    n_accepted = 0;
    n_tried    = 0;

    for iter = 1:max_iter_mc
        for i = 1:num_rods
            disp_vec = (rand - 0.5) * step_size * 2.0;

            e_curr    = overlap_energy_1D(i, centers(i), DD_current, centers, cell_head, cell_next, cell_size, ncell, L);
            trial_pos = mod(centers(i) + disp_vec, L);
            e_trial   = overlap_energy_1D(i, trial_pos, DD_current, centers, cell_head, cell_next, cell_size, ncell, L);

            n_tried = n_tried + 1;
            if e_trial <= e_curr
                % Remove from cell list
                ic = min(floor(centers(i) / cell_size) + 1, ncell);
                [cell_head, cell_next] = cell_remove_1D(i, ic, cell_head, cell_next);
                
                centers(i) = trial_pos;
                
                % Re-insert
                ic = min(floor(centers(i) / cell_size) + 1, ncell);
                cell_next(i) = cell_head(ic);
                cell_head(ic) = i;
                
                n_accepted = n_accepted + 1;
            end
        end
    end

    phi_current = (num_rods * DD_current) / L;
    n_overlaps  = count_overlaps_1D(num_rods, DD_current, centers, cell_head, cell_next, cell_size, ncell, L);
    fprintf('  DD/DD_target=%.4f  phi=%.4f  overlaps=%6d  acc%%=%.2f%%\n', ...
        DD_current/DD, phi_current, n_overlaps, 100*n_accepted/n_tried);

    if n_accepted/n_tried > 0.5
        step_size = step_size * 1.05;
    elseif n_accepted/n_tried < 0.2
        step_size = step_size * 0.95;
    end
    step_size = max(step_size, DD*1e-5);
    step_size = min(step_size, DD*0.5);
end

% =========================================================
% 4. Final hard-core relaxation
% =========================================================
fprintf('\n--- Phase 3: Final hard-core relaxation ---\n');

% IMPORTANT: Rebuild cell list with FINAL diameter
[cell_head, cell_next, cell_size, ncell] = build_cell_list_1D(DD, L, num_rods);
[cell_head, cell_next] = rebuild_cell_list_from_centers_1D(centers, num_rods, cell_size, ncell);

n_overlaps_resolved = count_overlaps_1D(num_rods, DD, centers, cell_head, cell_next, cell_size, ncell, L);
fprintf('  Initial overlaps at final DD: %d\n', n_overlaps_resolved);

step_size = DD * 0.01;

for phase = 1:2000
    n_overlaps_resolved = count_overlaps_1D(num_rods, DD, centers, cell_head, cell_next, cell_size, ncell, L);
    if n_overlaps_resolved == 0
        break;
    end

    n_accepted = 0;
    n_tried    = 0;

    for i = 1:num_rods
        e_curr = overlap_energy_1D(i, centers(i), DD, centers, cell_head, cell_next, cell_size, ncell, L);
        if e_curr < 1e-30
            continue;
        end

        best_pos = centers(i);
        best_e   = e_curr;
        for tr = 1:10
            disp_vec  = (rand - 0.5) * step_size * 2.0;
            trial_pos = mod(centers(i) + disp_vec, L);
            trial_e   = overlap_energy_1D(i, trial_pos, DD, centers, cell_head, cell_next, cell_size, ncell, L);
            n_tried   = n_tried + 1;
            if trial_e < best_e
                best_e   = trial_e;
                best_pos = trial_pos;
                n_accepted = n_accepted + 1;
            end
        end

        if best_pos ~= centers(i)
            ic = min(floor(centers(i)/cell_size)+1, ncell);
            [cell_head, cell_next] = cell_remove_1D(i, ic, cell_head, cell_next);
            centers(i) = best_pos;
            ic = min(floor(centers(i)/cell_size)+1, ncell);
            cell_next(i) = cell_head(ic);
            cell_head(ic) = i;
        end
    end

    if mod(phase, 100) == 0 || phase == 1
        fprintf('  Iter %4d  overlaps=%6d  step=%.4f\n', ...
            phase, n_overlaps_resolved, step_size/DD);
    end

    if n_tried > 0 && n_accepted/n_tried > 0.5
        step_size = step_size * 1.05;
    else
        step_size = step_size * 0.98;
    end
    step_size = max(step_size, DD*1e-6);
end

n_overlaps_resolved = count_overlaps_1D(num_rods, DD, centers, cell_head, cell_next, cell_size, ncell, L);
real_phi = (num_rods * DD) / L;
elapsed  = toc(t_start);

fprintf('\n=============================================\n');
fprintf(' Final Results\n');
fprintf('=============================================\n');
fprintf(' Rods              : %d\n', num_rods);
fprintf(' Final phi         : %g\n', real_phi);
fprintf(' Remaining overlaps: %d\n', n_overlaps_resolved);
fprintf(' CPU time (s)      : %.2f\n', elapsed);
if n_overlaps_resolved > 0
    fprintf(' WARNING: Packing not fully hard-core.\n');
end

end


%% ====================== 1D Helper Functions ======================

function [cell_head, cell_next, cell_size, ncell] = build_cell_list_1D(diam, L, num_rods)
    cell_size = max(diam, L/20);        % Prevent ncell=1 when diam is large
    ncell     = max(5, floor(L / cell_size));  % At least 5 cells when possible
    cell_head = zeros(ncell, 1, 'int32');
    cell_next = zeros(num_rods, 1, 'int32');
end

function [cell_head, cell_next] = rebuild_cell_list_from_centers_1D(centers, n, cell_size, ncell)
    cell_head = zeros(ncell, 1, 'int32');
    cell_next = zeros(n, 1, 'int32');
    for ii = 1:n
        icc = min(floor(centers(ii) / cell_size) + 1, ncell);
        cell_next(ii) = cell_head(icc);
        cell_head(icc) = ii;
    end
end

function [cell_head, cell_next] = cell_remove_1D(idx, icc, cell_head, cell_next)
    prev = 0;
    curr = cell_head(icc);
    while curr ~= 0 && curr ~= idx
        prev = curr;
        curr = cell_next(curr);
    end
    if curr == idx
        if prev == 0
            cell_head(icc) = cell_next(idx);
        else
            cell_next(prev) = cell_next(idx);
        end
        cell_next(idx) = 0;
    end
end

function energy = overlap_energy_1D(idx, pos, diam, centers, cell_head, cell_next, cell_size, ncell, L)
    energy = 0.0;
    icc = min(floor(pos / cell_size) + 1, ncell);
    for ddx = -1:1
        jcc = mod(icc + ddx - 1, ncell) + 1;
        nn  = cell_head(jcc);
        while nn ~= 0
            if nn ~= idx
                dd = abs(pos - centers(nn));
                if dd > L*0.5, dd = L - dd; end
                if dd < diam
                    energy = energy + (diam - dd)^2;
                end
            end
            nn = cell_next(nn);
        end
    end
end

function noverlap = count_overlaps_1D(n, diam, centers, cell_head, cell_next, cell_size, ncell, L)
    noverlap = 0;
    for ii = 1:n
        icc = min(floor(centers(ii) / cell_size) + 1, ncell);
        for ddx = -1:1
            jcc = mod(icc + ddx - 1, ncell) + 1;
            jj  = cell_head(jcc);
            while jj ~= 0
                if jj > ii
                    dd = abs(centers(ii) - centers(jj));
                    if dd > L*0.5, dd = L - dd; end
                    if dd < diam * (1.0 - 1e-10)
                        noverlap = noverlap + 1;
                    end
                end
                jj = cell_next(jj);
            end
        end
    end
end

