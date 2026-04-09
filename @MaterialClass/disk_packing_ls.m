function [centers, num_disks] = disk_packing_ls(L,DD,target_phi)
% disk_packing_ls.m
% 2D version of the RSA + Lubachevsky-Stillinger inspired grow-and-relax
% algorithm for disk packing in a square periodic box.

t_start = tic;

% =========================================================
% 1. Parameters
% =========================================================
escala     = 1.0;
% L          = [30e-3 50e-3] * escala;   % 2D box dimensions
% DD = DD*escala;
%DD         = 2.045e-2 * escala;      % target diameter
%target_phi = 0.45;

vol_disk   = pi * (DD/2)^2;          % area of one disk
vol_caixa  = L(1) * L(2);            % box area
num_disks  = ceil((vol_caixa * target_phi) / vol_disk);

fprintf('=============================================\n');
fprintf(' Disk Packing: RSA + LS Compression (2D)\n');
fprintf('=============================================\n');
fprintf(' Target phi        : %g\n', target_phi);
fprintf(' Target N disks    : %d\n', num_disks);
fprintf(' Box dimensions    : %g x %g\n', L(1), L(2));
fprintf(' Final diameter    : %g\n', DD);

centers = zeros(num_disks, 2);

rng('shuffle');

% =========================================================
% 2. Phase 1 — RSA at reduced diameter
% =========================================================
DD_current = DD * (0.42 / target_phi)^(1/2);   % 2D: exponent is 1/2
fprintf('\n--- Phase 1: RSA at reduced diameter ---\n');
fprintf('    RSA diameter    : %g\n', DD_current);
fprintf('    RSA phi target  : %g\n', ...
    (num_disks * pi*(DD_current/2)^2) / vol_caixa);

[cell_head, cell_next, cell_size, ncell] = build_cell_list(DD_current, L, num_disks);

max_att  = 500000000;
attempts = 0;
cnt      = 0;

while cnt < num_disks && attempts < max_att
    cand = rand(1,2) .* L;
    ic   = min(floor(cand ./ cell_size) + 1, ncell);

    overlapping = false;
    outerDone   = false;
    for dx = -1:1
        if outerDone, break; end
        for dy = -1:1
            jc = mod(ic + [dx, dy] - 1, ncell) + 1;
            ni = cell_head(jc(1), jc(2));
            while ni ~= 0
                d_vec = abs(cand - centers(ni,:));
                d_vec(d_vec > L*0.5) = d_vec(d_vec > L*0.5) - L(d_vec > L*0.5);
                dist2 = sum(d_vec.^2);
                if dist2 < DD_current^2
                    overlapping = true;
                    outerDone   = true;
                    break;
                end
                ni = cell_next(ni);
            end
            if outerDone, break; end
        end
    end

    if ~overlapping
        cnt = cnt + 1;
        centers(cnt,:) = cand;
        ic = min(floor(cand ./ cell_size) + 1, ncell);
        cell_next(cnt) = cell_head(ic(1), ic(2));
        cell_head(ic(1), ic(2)) = cnt;
        if mod(cnt, 1000) == 0
            fprintf('  RSA placed: %d / %d\n', cnt, num_disks);
        end
        attempts = 0;
    else
        attempts = attempts + 1;
    end
end

if cnt < num_disks
    error('ERROR: RSA could not place all disks at reduced diameter!');
end
fprintf('  RSA complete: %d disks placed.\n', cnt);

% =========================================================
% 3. Phase 2 — Grow diameters + MC displacement to relax
% =========================================================
fprintf('\n--- Phase 2: Grow & Relax (LS-MC) ---\n');

grow_rate   = 1.005;
max_iter_mc = 5000;
step_size   = DD * 0.05;

while DD_current < DD
    DD_current = min(DD_current * grow_rate, DD);

    [cell_head, cell_next, cell_size, ncell] = build_cell_list(DD_current, L, num_disks);
    [cell_head, cell_next] = rebuild_cell_list_from_centers(centers, num_disks, cell_size, ncell, cell_head, cell_next);

    n_accepted = 0;
    n_tried    = 0;

    for iter = 1:max_iter_mc
        for i = 1:num_disks
            disp_vec = (rand(1,2) - 0.5) * step_size * 2.0;

            e_curr    = overlap_energy(i, centers(i,:), DD_current, centers, cell_head, cell_next, cell_size, ncell, L);
            trial_pos = mod(centers(i,:) + disp_vec, L);
            e_trial   = overlap_energy(i, trial_pos, DD_current, centers, cell_head, cell_next, cell_size, ncell, L);

            n_tried = n_tried + 1;
            if e_trial <= e_curr
                ic = min(floor(centers(i,:) ./ cell_size) + 1, ncell);
                [cell_head, cell_next] = cell_remove(i, ic, cell_head, cell_next);
                centers(i,:) = trial_pos;
                ic = min(floor(centers(i,:) ./ cell_size) + 1, ncell);
                cell_next(i) = cell_head(ic(1), ic(2));
                cell_head(ic(1), ic(2)) = i;
                n_accepted = n_accepted + 1;
            end
        end
    end

    phi_current = (num_disks * pi*(DD_current/2)^2) / vol_caixa;
    n_overlaps  = count_overlaps(num_disks, DD_current, centers, cell_head, cell_next, cell_size, ncell, L);
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

[cell_head, cell_next] = rebuild_cell_list_from_centers(centers, num_disks, cell_size, ncell, cell_head, cell_next);
n_overlaps_resolved = count_overlaps(num_disks, DD, centers, cell_head, cell_next, cell_size, ncell, L);
fprintf('  Initial overlaps at final DD: %d\n', n_overlaps_resolved);

step_size = DD * 0.01;

for phase = 1:2000
    n_overlaps_resolved = count_overlaps(num_disks, DD, centers, cell_head, cell_next, cell_size, ncell, L);
    if n_overlaps_resolved == 0
        break;
    end

    n_accepted = 0;
    n_tried    = 0;

    for i = 1:num_disks
        e_curr = overlap_energy(i, centers(i,:), DD, centers, cell_head, cell_next, cell_size, ncell, L);
        if e_curr < 1e-30
            continue;
        end

        best_pos = centers(i,:);
        best_e   = e_curr;
        for tr = 1:10
            disp_vec  = (rand(1,2) - 0.5) * step_size * 2.0;
            trial_pos = mod(centers(i,:) + disp_vec, L);
            trial_e   = overlap_energy(i, trial_pos, DD, centers, cell_head, cell_next, cell_size, ncell, L);
            n_tried   = n_tried + 1;
            if trial_e < best_e
                best_e   = trial_e;
                best_pos = trial_pos;
                n_accepted = n_accepted + 1;
            end
        end

        if any(best_pos ~= centers(i,:))
            ic = min(floor(centers(i,:)./cell_size)+1, ncell);
            [cell_head, cell_next] = cell_remove(i, ic, cell_head, cell_next);
            centers(i,:) = best_pos;
            ic = min(floor(centers(i,:)./cell_size)+1, ncell);
            cell_next(i) = cell_head(ic(1), ic(2));
            cell_head(ic(1), ic(2)) = i;
        end
    end

    if mod(phase, 100) == 0
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

n_overlaps_resolved = count_overlaps(num_disks, DD, centers, cell_head, cell_next, cell_size, ncell, L);
real_phi = (num_disks * vol_disk) / vol_caixa;
elapsed  = toc(t_start);

fprintf('\n=============================================\n');
fprintf(' Final Results\n');
fprintf('=============================================\n');
fprintf(' Disks             : %d\n', num_disks);
fprintf(' Final phi         : %g\n', real_phi);
fprintf(' Remaining overlaps: %d\n', n_overlaps_resolved);
fprintf(' CPU time (s)      : %.2f\n', elapsed);
if n_overlaps_resolved > 0
    fprintf(' WARNING: Packing not fully hard-core. Consider more relaxation steps.\n');
end

%write_centers('disk_centers_45P_LSMETHOD.dat', centers, num_disks);

end % main function

% =========================================================
% Helper functions
% =========================================================

function [cell_head, cell_next, cell_size, ncell] = build_cell_list(diam, L, num_disks)
    cell_size = diam;
    ncell     = max(1, floor(L ./ cell_size));   % 1x2 vector
    cell_head = zeros(ncell(1), ncell(2), 'int32');
    cell_next = zeros(num_disks, 1, 'int32');
end

function [cell_head, cell_next] = rebuild_cell_list_from_centers(centers, n, cell_size, ncell, cell_head, cell_next)
    cell_head(:) = 0;
    cell_next(:) = 0;
    for ii = 1:n
        icc = min(floor(centers(ii,:) ./ cell_size) + 1, ncell);
        cell_next(ii) = cell_head(icc(1), icc(2));
        cell_head(icc(1), icc(2)) = ii;
    end
end

function [cell_head, cell_next] = cell_remove(idx, icc, cell_head, cell_next)
    prev = 0;
    curr = cell_head(icc(1), icc(2));
    while curr ~= 0 && curr ~= idx
        prev = curr;
        curr = cell_next(curr);
    end
    if curr == idx
        if prev == 0
            cell_head(icc(1), icc(2)) = cell_next(idx);
        else
            cell_next(prev) = cell_next(idx);
        end
        cell_next(idx) = 0;
    end
end

function energy = overlap_energy(idx, pos, diam, centers, cell_head, cell_next, cell_size, ncell, L)
    energy = 0.0;
    icc = min(floor(pos ./ cell_size) + 1, ncell);
    for ddx = -1:1
        for ddy = -1:1
            jcc = mod(icc + [ddx, ddy] - 1, ncell) + 1;
            nn  = cell_head(jcc(1), jcc(2));
            while nn ~= 0
                if nn ~= idx
                    dd2 = abs(pos - centers(nn,:));
                    dd2(dd2 > L*0.5) = dd2(dd2 > L*0.5) - L(dd2 > L*0.5);
                    dist2_val = sum(dd2.^2);
                    if dist2_val < diam^2
                        energy = energy + (diam - sqrt(dist2_val))^2;
                    end
                end
                nn = cell_next(nn);
            end
        end
    end
end

function noverlap = count_overlaps(n, diam, centers, cell_head, cell_next, cell_size, ncell, L)
    noverlap = 0;
    for ii = 1:n
        icc = min(floor(centers(ii,:) ./ cell_size) + 1, ncell);
        for ddx = -1:1
            for ddy = -1:1
                jcc = mod(icc + [ddx, ddy] - 1, ncell) + 1;
                jj  = cell_head(jcc(1), jcc(2));
                while jj ~= 0
                    if jj > ii
                        dd2 = abs(centers(ii,:) - centers(jj,:));
                        dd2(dd2 > L*0.5) = dd2(dd2 > L*0.5) - L(dd2 > L*0.5);
                        dist2_val = sum(dd2.^2);
                        if dist2_val < diam^2 * (1.0 - 1e-10)
                            noverlap = noverlap + 1;
                        end
                    end
                    jj = cell_next(jj);
                end
            end
        end
    end
end

function write_centers(filename, centers, n)
    fid = fopen(filename, 'w');
    fprintf(fid, '# x, y coordinates of disk centers\n');
    fprintf(fid, '# Number of disks: %d\n', n);
    for ii = 1:n
        fprintf(fid, '%23.15E %23.15E\n', centers(ii,1), centers(ii,2));
    end
    fclose(fid);
    fprintf('Centers saved to: %s\n', filename);
end

