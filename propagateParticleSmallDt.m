function P = propagateParticleSmallDt(mat,geometry,P,T)

% choice of time sub-step
T = T-mean(P.t);
dt = min(0.01*min(mat.meanFreeTime(:)),T); % with 0.01, we miss scattering events with a rate of 1e-5
Nt = ceil(T/dt);
dt = T/Nt;

% Pre-calculate Zoeppritz coefficients for elastic mode
% We do this outside the loop to avoid re-computing it at every sub-step
if ~mat.acoustics && isfield(geometry,'bnd')
    % Check if any boundary is cylindrical (dir==4) to skip
    hasCyl = false;
    for i2 = 1:length(geometry.bnd)
        if geometry.bnd(i2).dir == 4, hasCyl = true; break; end
    end
    
    if hasCyl
        % Fetch reflection coefficients
        [out_Z, ~] = MaterialClass.Zoeppritz(mat);
        
        % Create fast interpolants for the coefficients (Linear is usually sufficient/fast)
        % Note: Input angles for Zoeppritz are in degrees
        Rpp_Interp   = griddedInterpolant(out_Z.j1_deg, out_Z.E_Rpp,   'linear', 'nearest');
        Rsvsv_Interp = griddedInterpolant(out_Z.j1_deg, out_Z.E_Rsvsv, 'linear', 'nearest');
        % For solid-air (traction-free), energy conservation implies Rpp+Rpsv=1 (approx)
        % so we only strictly need Rpp to decide P vs S. Same for Rsvsv.
    end
end

% loop on time sub-iterations
for i1 = 1:Nt

    % polarization type of particles
    p = P.p;
    s = ~P.p;

    % propagate particles
    x1 = P.x;
    P.x(p,:) = P.x(p,:) + (mat.vp.*dt).*P.dir(p,:);
    P.x(s,:) = P.x(s,:) + (mat.vs.*dt).*P.dir(s,:);

    % bounce on boundaries
    if isfield(geometry,'bnd')

        for i2 = 1:length(geometry.bnd)

            bnd = geometry.bnd(i2);

            % boundaries along cartesian coordinates
            if bnd.dir<4
                % index of particles outside the boundary
                ind = (P.x(:,bnd.dir)-bnd.val).*(x1(:,bnd.dir)-bnd.val)<0;
                
                if ~any(ind); continue; end
                
                if mat.acoustics
                    % Acoustic: specular reflection
                    P.x(ind,bnd.dir) = (2*bnd.val)-P.x(ind,bnd.dir);
                    P.dir(ind,bnd.dir) = -P.dir(ind,bnd.dir);
                else
                    % Elastic with plane boundary
                    n_plane = zeros(1, P.d); 
                    n_plane(bnd.dir) = 1; % unit normal to the interface

                    % Decide P/S and SV/SH behavior (reflection coeffs)
                    % MaterialInterface handles the physics/probabilities
                    o = mat.MaterialInterface(P, n_plane, ind);

                    % Reflect positions (only for plane boundaries)
                    P.x(ind,bnd.dir) = (2*bnd.val)-P.x(ind,bnd.dir);

                    % Reflection & Mode conversion
                    
                    % -- P Incident --
                    % P -> P
                    if any(o.p2p)
                        idx = o.pparticle;   % P particles hitting boundary
                        idx(idx) = o.p2p;    % refine to P->P
                        P.dir(idx, bnd.dir)  = -P.dir(idx,  bnd.dir);
                        P.perp(idx, bnd.dir) = -P.perp(idx, bnd.dir);
                    end
                    
                    % P -> S 
                    if any(o.p2sv)
                        idx = o.pparticle;   % P particles hitting boundary
                        idx(idx) = o.p2sv;   % refine to P->S
                        
                        % Update polarization flag to S (false)
                        P.p(idx) = false;
                        
                        % output-to-input velocity ratio
                        velocityRatio = mat.vs / mat.vp; 
                        
                        % Tangent vector conservation
                        d_old = P.dir(idx, :);
                        d_n = d_old(:, bnd.dir); % Normal component
                        d_t = d_old; 
                        d_t(:, bnd.dir) = 0;     % Tangent component

                        d_new_t = velocityRatio * d_t; % Snell's Law
                        
                        % New normal Component (pointing away from boundary)
                        sin2 = sum(d_new_t.^2, 2);
                        cos_val = sqrt(max(0, 1 - sin2));
                        
                        d_new = d_new_t;
                        d_new(:, bnd.dir) = -sign(d_n) .* cos_val;
                        
                        % Normalize
                        P.dir(idx, :) = d_new ./ vecnorm(d_new, 2, 2);
                        
                        % Update Perp : 
                        % SV polarization is in plane of incidence that 
                        % contains n and d. n is aligned with axis bnd.dir.
                        n_vec = zeros(sum(idx), P.d);
                        n_vec(:, bnd.dir) = 1;
                        axis_rot = cross(P.dir(idx, :), n_vec, 2);
                        axis_rot = axis_rot ./ vecnorm(axis_rot, 2, 2);
                        P.perp(idx, :) = cross(axis_rot, P.dir(idx, :), 2);
                    end
                    
                    % -- S Incident --
                    % SH -> SH
                    if any(o.sh2sh)
                         idx = o.sparticle;       % S particles hitting boundary
                         idx(idx) = o.shparticle; % refine to SH particles
                         P.dir(idx, bnd.dir)  = -P.dir(idx,  bnd.dir);
                         P.perp(idx, bnd.dir) = -P.perp(idx, bnd.dir);
                    end

                    % SV -> SV
                    if any(o.sv2sv)
                         % Only take those that are SV AND reflected as SV
                         isSV = false(P.N,1);
                         isSV(o.sparticle) = o.svparticle; % Mask of all SV
                         
                         % intersect with decision
                         isSVtoP = false(P.N,1);
                         isSVtoP(isSV) = o.sv2sv;
                         
                         P.dir(isSVtoP,  bnd.dir) = -P.dir(isSVtoP, bnd.dir);
                         P.perp(isSVtoP, bnd.dir) = -P.perp(isSVtoP, bnd.dir);
                    end
                    
                    % SV -> P
                    if any(o.sv2p)
                         isSV = false(P.N,1);
                         % TRUE for SV particles at the boundary
                         isSV(o.sparticle) = o.svparticle;

                         isSVtoP = false(P.N,1);
                         % TRUE for SV particles that are converting to P at the boundary
                         isSVtoP(isSV) = o.sv2p;
                         
                         if any(isSVtoP)
                             % Update polarization flag to P (true)
                             P.p(isSVtoP) = true;
                             
                             velocityRatio = mat.vp / mat.vs;
                             
                             d_old = P.dir(isSVtoP, :);
                             d_n = d_old(:, bnd.dir);
                             d_t = d_old; d_t(:, bnd.dir) = 0;
                             
                             d_new_t = velocityRatio * d_t; % Snell's Law
                             sin2 = sum(d_new_t.^2, 2);
                             
                             % Handle Critical Angle (Total Internal Reflection)
                             % If sin2 > 1, no transmission to P (evanescent)
                             % Reflect as SV instead
                             valid = sin2 <= 1;
                             
                             if any(valid)
                                 % Update valid conversions
                                 sub_idx = find(isSVtoP);
                                 valid_idx = sub_idx(valid);
                                 
                                 cos_val = sqrt(1 - sin2(valid));
                                 d_new = d_new_t(valid, :);
                                 d_new(:, bnd.dir) = -sign(d_n(valid)) .* cos_val;
                                 
                                 P.dir(valid_idx, :) = d_new ./ vecnorm(d_new, 2, 2);
                             end
                             
                             if any(~valid)
                                 % Revert invalid conversions (TIR -> SV Reflection)
                                 sub_idx = find(isSVtoP);
                                 invalid_idx = sub_idx(~valid);
                                 P.p(invalid_idx) = false; % Revert to S
                                 P.dir(invalid_idx, bnd.dir) = -P.dir(invalid_idx, bnd.dir);
                             end
                         end
                    end
                end

            % boundary along cylindrical radius
            elseif bnd.dir==4
                % distance to z axis
                radial_coords = P.x(:, 1:2);
                r = sqrt(sum(radial_coords.^2, 2));

                % which particles have moved outside the cylinder radius
                overshoot_dist = r - bnd.val;
                ind = overshoot_dist > 0;

                if any(ind)
                    pos_outside = P.x(ind, :);
                    r_outside = r(ind);
            
                    % unit normal vector (in x-y plane) at impact point
                    % Points OUTWARDS from cylinder center
                    n = zeros(size(pos_outside));
                    n(:, 1:2) = pos_outside(:, 1:2) ./ r_outside;

                    % update position (Geometric Mirroring - Same for Acoustic/Elastic)
                    P.x(ind, :) = pos_outside - 2 * overshoot_dist(ind) .* n;

                    if mat.acoustics
                        % Acoustic: Specular Reflection
                        % update direction
                        dir_outside = P.dir(ind, :);
                        dot_prod = sum(dir_outside .* n, 2);
                        P.dir(ind,:) = dir_outside - 2 * dot_prod .* n;
                    else
                        % Elastic: Curved Boundary Reflection
                        % Note: We cannot simply use mat.MaterialInterface because 
                        % n varies per particle. We must implement the logic inline.
                        
                        % Indices of particles hitting the boundary
                        idx_hit = find(ind);
                        
                        % --- 1. Identify P vs S incident ---
                        isP = P.p(ind);
                        isS = ~P.p(ind);
                        
                        % --- 2. Handle P Incident Particles ---
                        if any(isP)
                            % Global indices
                            glob_p = idx_hit(isP);
                            
                            % Vectors for P particles
                            n_p = n(isP, :);
                            d_p = P.dir(glob_p, :);
                            p_p = P.perp(glob_p, :); % Perp vector
                            
                            % Calculate Incident Angle (dot product)
                            % d points out, n points out -> dot > 0
                            cos_theta = sum(d_p .* n_p, 2);
                            theta_deg = acosd(min(1, max(-1, cos_theta)));
                            
                            % Get Rpp coefficient
                            Rpp_val = Rpp_Interp(theta_deg);
                            
                            % Monte Carlo: P->P or P->S?
                            rand_val = rand(length(glob_p), 1);
                            mask_p2p = rand_val <= Rpp_val;
                            mask_p2s = ~mask_p2p;
                            
                            % A. Specular Reflection (P->P)
                            if any(mask_p2p)
                                g_idx = glob_p(mask_p2p);
                                n_sub = n_p(mask_p2p, :);
                                d_sub = d_p(mask_p2p, :);
                                p_sub = p_p(mask_p2p, :);
                                
                                % Reflect Direction: d - 2(d.n)n
                                dt_n = sum(d_sub .* n_sub, 2);
                                P.dir(g_idx, :) = d_sub - 2 * dt_n .* n_sub;
                                
                                % Reflect Perp (keep orthogonal)
                                pt_n = sum(p_sub .* n_sub, 2);
                                P.perp(g_idx, :) = p_sub - 2 * pt_n .* n_sub;
                            end
                            
                            % B. Mode Conversion (P->S)
                            if any(mask_p2s)
                                g_idx = glob_p(mask_p2s);
                                n_sub = n_p(mask_p2s, :);
                                d_sub = d_p(mask_p2s, :);
                                
                                % Change type
                                P.p(g_idx) = false;
                                
                                % output-to-input velocity ratio
                                eta_val = mat.vs / mat.vp;
                                
                                % Decompose: d = d_n * n + d_t
                                dot_dn = sum(d_sub .* n_sub, 2);
                                d_n_vec = dot_dn .* n_sub;
                                d_t_vec = d_sub - d_n_vec;
                                
                                % New Tangent: d_new_t = eta * d_t
                                d_new_t = eta_val * d_t_vec;
                                
                                % New Normal Magnitude
                                sin2 = sum(d_new_t.^2, 2);
                                cos_val = sqrt(max(0, 1 - sin2));
                                
                                % New Direction: d_new = d_new_t - cos * n
                                % (Reflected wave goes IN, so against n)
                                d_new = d_new_t - cos_val .* n_sub;
                                d_new = d_new ./ vecnorm(d_new, 2, 2);
                                P.dir(g_idx, :) = d_new;
                                
                                % Update Perp (SV polarization in plane of incidence)
                                % Plane defined by (d_new, n). 
                                % Axis = d_new x n. SV = axis x d_new.
                                axis_vec = cross(d_new, n_sub, 2);
                                axis_len = vecnorm(axis_vec, 2, 2);
                                % Handle normal incidence singularity
                                mask_sing = axis_len < 1e-6;
                                if any(mask_sing)
                                    % Arbitrary perp if normal incidence
                                    axis_vec(mask_sing, :) = cross(d_new(mask_sing,:), repmat([1 0 0], sum(mask_sing), 1), 2);
                                end
                                axis_vec = axis_vec ./ vecnorm(axis_vec, 2, 2);
                                P.perp(g_idx, :) = cross(axis_vec, d_new, 2);
                            end
                        end
                        
                        % --- 3. Handle S Incident Particles ---
                        if any(isS)
                            glob_s = idx_hit(isS);
                            n_s = n(isS, :);
                            d_s = P.dir(glob_s, :);
                            p_s = P.perp(glob_s, :);
                            
                            % Decompose into SV and SH
                            % Plane of incidence normal (SH direction): axis = d x n
                            axis_inc = cross(d_s, n_s, 2);
                            axis_inc = axis_inc ./ vecnorm(axis_inc, 2, 2);
                            
                            % Project polarization onto SH axis
                            sh_comp = dot(p_s, axis_inc, 2);
                            
                            % Probability of being SH (Energy fraction)
                            prob_sh = sh_comp.^2; 
                            
                            rand_type = rand(length(glob_s), 1);
                            mask_sh = rand_type < prob_sh;
                            mask_sv = ~mask_sh;
                            
                            % A. SH Particles (Specular Reflection)
                            if any(mask_sh)
                                g_idx = glob_s(mask_sh);
                                n_sub = n_s(mask_sh, :);
                                d_sub = d_s(mask_sh, :);
                                p_sub = p_s(mask_sh, :);
                                
                                % Reflect d
                                dt_n = sum(d_sub .* n_sub, 2);
                                P.dir(g_idx, :) = d_sub - 2 * dt_n .* n_sub;
                                
                                % Reflect p
                                pt_n = sum(p_sub .* n_sub, 2);
                                P.perp(g_idx, :) = p_sub - 2 * pt_n .* n_sub;
                            end
                            
                            % B. SV Particles (Specular or Conversion)
                            if any(mask_sv)
                                g_idx = glob_s(mask_sv);
                                n_sub = n_s(mask_sv, :);
                                d_sub = d_s(mask_sv, :);
                                
                                % Incident Angle
                                cos_theta = sum(d_sub .* n_sub, 2);
                                theta_deg = acosd(min(1, max(-1, cos_theta)));
                                
                                % Get Rsvsv
                                Rsvsv_val = Rsvsv_Interp(theta_deg);
                                
                                rand_conv = rand(length(g_idx), 1);
                                mask_sv2sv = rand_conv <= Rsvsv_val;
                                mask_sv2p = ~mask_sv2sv;
                                
                                % SV->SV (Specular)
                                if any(mask_sv2sv)
                                    g_fin = g_idx(mask_sv2sv);
                                    n_fin = n_sub(mask_sv2sv, :);
                                    d_fin = d_sub(mask_sv2sv, :);
                                    p_fin = P.perp(g_fin, :); % Original perp
                                    
                                    dt_n = sum(d_fin .* n_fin, 2);
                                    P.dir(g_fin, :) = d_fin - 2 * dt_n .* n_fin;
                                    
                                    pt_n = sum(p_fin .* n_fin, 2);
                                    P.perp(g_fin, :) = p_fin - 2 * pt_n .* n_fin;
                                end
                                
                                % SV->P (Conversion)
                                if any(mask_sv2p)
                                    g_fin = g_idx(mask_sv2p);
                                    n_fin = n_sub(mask_sv2p, :);
                                    d_fin = d_sub(mask_sv2p, :);
                                    
                                    P.p(g_fin) = true;
                                    
                                    eta_val = mat.vp / mat.vs;
                                    
                                    dot_dn = sum(d_fin .* n_fin, 2);
                                    d_n_vec = dot_dn .* n_fin;
                                    d_t_vec = d_fin - d_n_vec;
                                    
                                    d_new_t = eta_val * d_t_vec;
                                    sin2 = sum(d_new_t.^2, 2);
                                    
                                    % Check Critical Angle (Total Internal Reflection)
                                    valid = sin2 <= 1;
                                    
                                    if any(valid)
                                        g_val = g_fin(valid);
                                        n_val = n_fin(valid, :);
                                        d_t_val = d_new_t(valid, :);
                                        cos_val = sqrt(1 - sin2(valid));
                                        
                                        d_new = d_t_val - cos_val .* n_val;
                                        P.dir(g_val, :) = d_new ./ vecnorm(d_new, 2, 2);
                                        
                                        % P polarization (usually ignored, but good to set)
                                        % Just keep perp orthogonal or rotate it
                                    end
                                    
                                    if any(~valid)
                                        % TIR fallback -> SV reflection
                                        g_inv = g_fin(~valid);
                                        n_inv = n_fin(~valid, :);
                                        d_inv = P.dir(g_inv, :); % Still old dir
                                        p_inv = P.perp(g_inv, :);
                                        
                                        P.p(g_inv) = false; % Revert type
                                        
                                        dt_n = sum(d_inv .* n_inv, 2);
                                        P.dir(g_inv, :) = d_inv - 2 * dt_n .* n_inv;
                                        
                                        pt_n = sum(p_inv .* n_inv, 2);
                                        P.perp(g_inv, :) = p_inv - 2 * pt_n .* n_inv;
                                    end
                                end
                            end
                        end
                        % End Elastic Logic
                    end
                end
            end
        end
    end

    % choose particles that are scattered at the end of time step
    if isscalar(mat.meanFreeTime)
        mat.meanFreeTime = [mat.meanFreeTime; mat.meanFreeTime];
    end

    lambda = zeros(P.N,1);
    lambda(p) = dt / mat.meanFreeTime(1); 
    lambda(s) = dt / mat.meanFreeTime(2);
    
    Nscat = poissrnd(lambda);
    scat  = Nscat > 0;
    Nscat = sum(scat);

    % draw scattering angle around direction of propagation
    if P.d==2
        phi = randi([0 1],Nscat,1,'logical');
        P.perp(scat,:) = ((-1).^phi).*P.perp(scat,:);
    elseif P.d==3
        phi = (2*pi)*rand(Nscat,1);
        P.perp(scat,:) = cos(phi).*P.perp(scat,:) ...
            + sin(phi).*cross(P.dir(scat,:),P.perp(scat,:));
    end

    % draw scattering angle away from direction of propagation
    if mat.acoustics
        theta = zeros(P.N,1);
        theta(scat) = mat.invcdf(rand(Nscat,1));

    else
        % change polarization
        p = scat & p;
        P.p( p & (rand(P.N,1)>mat.P2P) ) = false;
        s = scat & s;
        P.p( s & (rand(P.N,1)>mat.S2S) ) = true;

        % draw scattering angle away from direction of propagation
        theta = zeros(P.N,1);
        theta( scat ) = rand(Nscat,1);
        theta(  P.p & p ) = mat.invcdf{1,1}( theta(  P.p & p ) );
        theta( ~P.p & p ) = mat.invcdf{1,2}( theta( ~P.p & p ) );
        theta(  P.p & s ) = mat.invcdf{2,1}( theta(  P.p & s ) );
        theta( ~P.p & s ) = mat.invcdf{2,2}( theta( ~P.p & s ) );
    end

    % scattering away from direction of propagation
    dir = P.dir;
    P.dir =  cos(theta).*dir + sin(theta).*P.perp;
    P.dir = P.dir./sqrt(sum(P.dir.^2,2));
    P.perp = -sin(theta).*dir + cos(theta).*P.perp;
    P.perp = P.perp./sqrt(sum(P.perp.^2,2));

    % end of loop on time sub-iterations
end

% update particle times
P.t = P.t+T;

end