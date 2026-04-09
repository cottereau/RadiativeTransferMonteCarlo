function [centers, num_spheres] = sphere_packing_ls(L,DD,target_phi)
% RSA placement followed by Lubachevsky-Stillinger inspired grow-and-relax
% compression to reach a target volume fraction (phi).

t_start = tic;

% =========================================================
% 1. Parameters
% =========================================================
pi_val     = pi;
escala     = 1.0;
% L          = [1.0, 1.0, 1.0] * escala;
% DD         = 2.045e-2 * escala;   % target diameter
% target_phi = 0.25;

vol_esfera  = (4/3) * pi_val * (DD/2)^3;
vol_caixa   = L(1) * L(2) * L(3);
num_spheres = ceil((vol_caixa * target_phi) / vol_esfera);

fprintf('=============================================\n');
fprintf(' Sphere Packing: RSA + LS Compression\n');
fprintf('=============================================\n');
fprintf(' Target phi        : %g\n', target_phi);
fprintf(' Target N spheres  : %d\n', num_spheres);
fprintf(' Box dimensions    : %g %g %g\n', L(1), L(2), L(3));
fprintf(' Final diameter    : %g\n', DD);

centers = zeros(num_spheres, 3);

% Seed random number generator with clock-based seed
rng('shuffle');

% =========================================================
% 2. Phase 1 — RSA at reduced diameter
%    Start with DD_rsa such that phi_rsa ~ 0.37 (easily reachable)
% =========================================================
DD_current = DD * (0.3 / target_phi)^(1/3);
fprintf('\n--- Phase 1: RSA at reduced diameter ---\n');
fprintf('    RSA diameter    : %g\n', DD_current);
fprintf('    RSA phi target  : %g\n', ...
    (num_spheres * (4/3)*pi_val*(DD_current/2)^3) / vol_caixa);

[cell_head, cell_next, cell_size, ncell] = build_cell_list(DD_current, L, num_spheres);

max_att  = 500000000;
attempts = 0;
cnt      = 0;

while cnt < num_spheres && attempts < max_att
    cand = rand(1,3) .* L;
    ic   = min(floor(cand ./ cell_size) + 1, ncell);

    overlapping = false;
    outerDone = false;
    for dx = -1:1
        if outerDone, break; end
        for dy = -1:1
            if outerDone, break; end
            for dz = -1:1
                jc = mod(ic + [dx, dy, dz] - 1, ncell) + 1;
                ni = cell_head(jc(1), jc(2), jc(3));
                while ni ~= 0
                    d_vec = abs(cand - centers(ni,:));
                    d_vec(d_vec > L*0.5) = d_vec(d_vec > L*0.5) - L(d_vec > L*0.5);
                    dist2 = sum(d_vec.^2);
                    if dist2 < DD_current^2
                        overlapping = true;
                        outerDone = true;
                        break;
                    end
                    ni = cell_next(ni);
                end
                if outerDone, break; end
            end
        end
    end

    if ~overlapping
        cnt = cnt + 1;
        centers(cnt,:) = cand;
        ic = min(floor(cand ./ cell_size) + 1, ncell);
        cell_next(cnt) = cell_head(ic(1), ic(2), ic(3));
        cell_head(ic(1), ic(2), ic(3)) = cnt;
        if mod(cnt, 1000) == 0
            fprintf('  RSA placed: %d / %d\n', cnt, num_spheres);
        end
        attempts = 0;
    else
        attempts = attempts + 1;
    end
end

if cnt < num_spheres
    error('ERROR: RSA could not place all spheres at reduced diameter!');
end
fprintf('  RSA complete: %d spheres placed.\n', cnt);

% =========================================================
% 3. Phase 2 — Grow diameters + MC displacement to relax
%    Lubachevsky-Stillinger inspired
% =========================================================
fprintf('\n--- Phase 2: Grow & Relax (LS-MC) ---\n');

grow_rate   = 1.005;
max_iter_mc = 5000;
step_size   = DD * 0.05;

while DD_current < DD
    DD_current = min(DD_current * grow_rate, DD);

    [cell_head, cell_next, cell_size, ncell] = build_cell_list(DD_current, L, num_spheres);
    [cell_head, cell_next] = rebuild_cell_list_from_centers(centers, num_spheres, cell_size, ncell, cell_head, cell_next);

    n_accepted = 0;
    n_tried    = 0;

    for iter = 1:max_iter_mc
        for i = 1:num_spheres
            disp_vec = (rand(1,3) - 0.5) * step_size * 2.0;

            e_curr = overlap_energy(i, centers(i,:), DD_current, centers, cell_head, cell_next, cell_size, ncell, L);

            trial_pos = centers(i,:) + disp_vec;
            trial_pos = mod(trial_pos, L);
            e_trial   = overlap_energy(i, trial_pos, DD_current, centers, cell_head, cell_next, cell_size, ncell, L);

            n_tried = n_tried + 1;
            if e_trial <= e_curr
                ic = min(floor(centers(i,:) ./ cell_size) + 1, ncell);
                [cell_head, cell_next] = cell_remove(i, ic, cell_head, cell_next);
                centers(i,:) = trial_pos;
                ic = min(floor(centers(i,:) ./ cell_size) + 1, ncell);
                cell_next(i) = cell_head(ic(1), ic(2), ic(3));
                cell_head(ic(1), ic(2), ic(3)) = i;
                n_accepted = n_accepted + 1;
            end
        end
    end

    phi_current = (num_spheres * (4/3)*pi_val*(DD_current/2)^3) / vol_caixa;
    n_overlaps  = count_overlaps(num_spheres, DD_current, centers, cell_head, cell_next, cell_size, ncell, L);
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
% 4. Final hard-core relaxation (zero overlaps required)
% =========================================================
fprintf('\n--- Phase 3: Final hard-core relaxation ---\n');

[cell_head, cell_next] = rebuild_cell_list_from_centers(centers, num_spheres, cell_size, ncell, cell_head, cell_next);
n_overlaps_resolved = count_overlaps(num_spheres, DD, centers, cell_head, cell_next, cell_size, ncell, L);
fprintf('  Initial overlaps at final DD: %d\n', n_overlaps_resolved);

step_size = DD * 0.01;

for phase = 1:2000
    n_overlaps_resolved = count_overlaps(num_spheres, DD, centers, cell_head, cell_next, cell_size, ncell, L);
    if n_overlaps_resolved == 0
        break;
    end

    n_accepted = 0;
    n_tried    = 0;

    for i = 1:num_spheres
        e_curr = overlap_energy(i, centers(i,:), DD, centers, cell_head, cell_next, cell_size, ncell, L);
        if e_curr < 1e-30
            continue;
        end

        best_pos = centers(i,:);
        best_e   = e_curr;
        for tr = 1:10
            disp_vec  = (rand(1,3) - 0.5) * step_size * 2.0;
            trial_pos = centers(i,:) + disp_vec;
            trial_pos = mod(trial_pos, L);
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
            cell_next(i) = cell_head(ic(1), ic(2), ic(3));
            cell_head(ic(1), ic(2), ic(3)) = i;
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

n_overlaps_resolved = count_overlaps(num_spheres, DD, centers, cell_head, cell_next, cell_size, ncell, L);
real_phi = (num_spheres * vol_esfera) / vol_caixa;
elapsed  = toc(t_start);

fprintf('\n=============================================\n');
fprintf(' Final Results\n');
fprintf('=============================================\n');
fprintf(' Spheres           : %d\n', num_spheres);
fprintf(' Final phi         : %g\n', real_phi);
fprintf(' Remaining overlaps: %d\n', n_overlaps_resolved);
fprintf(' CPU time (s)      : %.2f\n', elapsed);
if n_overlaps_resolved > 0
    fprintf(' WARNING: Packing not fully hard-core. Consider more relaxation steps.\n');
end

%write_centers('sphere_centers_45P_LSMETHOD_Herve.dat', centers, num_spheres);

end % main function

% =========================================================
% Helper functions
% =========================================================

function [cell_head, cell_next, cell_size, ncell] = build_cell_list(diam, L, num_spheres)
    cell_size = diam;
    ncell     = max(1, floor(L ./ cell_size));
    cell_head = zeros(ncell(1), ncell(2), ncell(3), 'int32');
    cell_next = zeros(num_spheres, 1, 'int32');
end

function [cell_head, cell_next] = rebuild_cell_list_from_centers(centers, n, cell_size, ncell, cell_head, cell_next)
    cell_head(:) = 0;
    cell_next(:) = 0;
    for ii = 1:n
        icc = min(floor(centers(ii,:) ./ cell_size) + 1, ncell);
        cell_next(ii) = cell_head(icc(1), icc(2), icc(3));
        cell_head(icc(1), icc(2), icc(3)) = ii;
    end
end

function [cell_head, cell_next] = cell_remove(idx, icc, cell_head, cell_next)
    prev = 0;
    curr = cell_head(icc(1), icc(2), icc(3));
    while curr ~= 0 && curr ~= idx
        prev = curr;
        curr = cell_next(curr);
    end
    if curr == idx
        if prev == 0
            cell_head(icc(1), icc(2), icc(3)) = cell_next(idx);
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
            for ddz = -1:1
                jcc = mod(icc + [ddx, ddy, ddz] - 1, ncell) + 1;
                nn  = cell_head(jcc(1), jcc(2), jcc(3));
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
end

function noverlap = count_overlaps(n, diam, centers, cell_head, cell_next, cell_size, ncell, L)
    noverlap = 0;
    for ii = 1:n
        icc = min(floor(centers(ii,:) ./ cell_size) + 1, ncell);
        for ddx = -1:1
            for ddy = -1:1
                for ddz = -1:1
                    jcc = mod(icc + [ddx, ddy, ddz] - 1, ncell) + 1;
                    jj  = cell_head(jcc(1), jcc(2), jcc(3));
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
end

function write_centers(filename, centers, n)
    fid = fopen(filename, 'w');
    fprintf(fid, '# x, y, z coordinates of sphere centers\n');
    fprintf(fid, '# Number of spheres: %d\n', n);
    for ii = 1:n
        fprintf(fid, '%23.15E %23.15E %23.15E\n', centers(ii,1), centers(ii,2), centers(ii,3));
    end
    fclose(fid);
    fprintf('Centers saved to: %s\n', filename);
end