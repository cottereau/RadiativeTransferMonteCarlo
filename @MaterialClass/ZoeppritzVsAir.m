function [out, angles] = ZoeppritzVsAir(mat)
    %% ZoeppritzVsAir
    % function to calculate reflection and transmission
    % coefficients for P, SV and SH waves
    %
    % Syntax:
    %   MaterialClass.ZoeppritzVsAir( mat );
    %
    % Inputs:
    %  mat: scalar elastic MaterialClass object. The exterior
    %       medium is currently hard-coded as air.
    %
    % Output: The method return the coefficients of transmission
    % and reflection or only reflection
    %
    % IMPORTANT -- this is NOT a free surface.  For a scalar mat the
    % second medium is hard-coded AIR (vp = 343 m/s, vs = 0,
    % rho = 1 kg/m^3), so a small but nonzero amount of energy is
    % genuinely transmitted out of the solid and the reflection
    % coefficients do NOT sum to 1 by design:
    %
    %   1 - (E_Rpp + E_Rpsv) = E_Tpp, angle for angle
    %
    % At normal incidence that equals the elementary acoustic
    % transmission coefficient 4*Z1*Z2/(Z1+Z2)^2 with Z1 = rho*vp and
    % Z2 = 1*343.  Measured for vp = 2500, vs = vp/sqrt(3),
    % rho = 2000 (Examples/BoundaryReflectionEnergyTests.m):
    %
    %   theta    1-(E_Rpp+E_Rpsv)   E_Tpp        ratio
    %     0 deg    2.7436e-04       2.7436e-04   1.0000
    %    15 deg    2.6156e-04       2.6156e-04   1.0000
    %    30 deg    2.2677e-04       2.2677e-04   1.0000
    %    45 deg    1.8044e-04       1.8044e-04   1.0000
    %    60 deg    1.3816e-04       1.3816e-04   1.0000
    %
    % E_Rsh sums to 1 exactly, because a fluid carries no SH.
    %
    % Callers that treat the boundary as perfectly reflecting must
    % decide where that ~2.7e-04 goes.  The RTE's free-surface path
    % gives it to the converted (SV) channel by sampling against
    % 1 - E_Rpp rather than E_Rpsv.  Small and deliberate, but it is
    % why a reflecting-boundary energy budget closes at 1 while these
    % raw coefficients do not.

    if isscalar(mat)
        % assuming solid-fluid interface
        % the fluid is hard coded as air
        vpAir  = 343; %m/s
        vsAir  = 0;   %m/s
        rhoAir = 1;  %kg/m^3

        % angles
        j1_deg = linspace(0,90,181);

        if(mat.acoustics)
            error("This dont work for acoustic material")
        end

        vp1  = mat(1).vp;
        vs1  = mat(1).vs;
        rho1 = mat(1).rho;

        vp2  = vpAir;
        vs2  = vsAir;
        rho2 = rhoAir;

        if vs2 == 0
            out = MaterialClass.ZoeppritzFluid(j1_deg,vp1,vs1,rho1,vp2,vs2,rho2);
        else
            out = MaterialClass.ZoeppritzSolid(j1_deg,vp1,vs1,rho1,vp2,vs2,rho2);
        end

        out.j1_deg = j1_deg;

        % angles P
        angles.Rpp  = @(theta1) asind((sind(theta1) / vp1) * vp1);
        angles.Rpsv = @(theta1) asind((sind(theta1) / vp1) * vs1);
        angles.Tpp  = @(theta1) asind((sind(theta1) / vp1) * vp2);
        angles.Tpsv = @(theta1) asind((sind(theta1) / vp1) * vs2);

        % angles SV
        angles.Rsvsv = @(theta1) asind(sind(theta1) / vs1 * vs1);
        angles.Rsvp  = @(theta1) asind(sind(theta1) / vs1 * vp1);
        angles.Tsvsv = @(theta1) asind(sind(theta1) / vs1 * vs2);
        angles.Tsvp  = @(theta1) asind(sind(theta1) / vs1 * vp2);

        % Incident SH-wave
        angles.Rsh = @(theta1) asind(sind(theta1) / vs1 * vs1);
        angles.Tsh = @(theta1) asind(sind(theta1) / vs1 * vs2);

    else
        error('MaterialClass:MaterialInterfacesNotImplemented', ...
            ['ZoeppritzVsAir accepts only one material and ', ...
            'models its interface with air. Interfaces between ', ...
            'multiple user-defined materials are not implemented.']);
    end
end
