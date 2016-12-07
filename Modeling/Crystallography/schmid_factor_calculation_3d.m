function [TPSCHMID, TPSCHMIDALL, TVSCHMID, TVSCHMIDALL, PSCHMID, PSCHMIDALL] = schmid_factor_calculation_3d(texture_ip, schmid_direction_ip, loading_type_ip, U_ip)
%% Calculate Schmid factor for slip, phase transformation
% Inputs:
%   texture_ip: (Nx4) crystal orientation in terms of quaternions
%   schmid_direction_ip: (int) direction for Schmid factor calculation.
%       This is the loading axis direction. x = 1, y = 2, z = 3
%   loading_type_ip: (int) +1 for tension, -1 for compression
%   U_ip: 
%
% Schmid factor calculation
schmid_direction = schmid_direction_ip; % x = 1, y = 2, z = 3
loading_type = loading_type_ip;         % -1 = compression, +1 = tension
%
U = U_ip;
num_variants = numel(U);
%
E = cell(num_variants,1);
for i=1:num_variants
    E{i} = 0.5*(transpose(U{i})*U{i} - eye(3));
end
% Texture as quaternion array
TEXTURE = texture_ip;
%
% In TEXTURE each row is an Euler triplet.
TEXSIZE = size(TEXTURE, 1);
PSCHMID = zeros(TEXSIZE, 1);
TVSCHMID = zeros(TEXSIZE, 1);
TPSCHMID = zeros(TEXSIZE, 1);
%
[cb0t, cm0t] = habit_calculation_3d(U);
% Loop over all Euler triplets.
for itex=1:TEXSIZE
    quat = TEXTURE(itex, :);
    tlg=quat2g(quat);
    % Schmid factor for plasticity assuming cubic
    % 12 slip systems. Right now 6 gives [110] systems.
    nslip = 12;
    % Get b, m in crystal coordinate system
    [cb0,cm0] = slipsystem_p(nslip);
    % Transform to sample
    b0 = tlg*cb0;
    m0 = tlg*cm0;
    s0alpha = zeros(3,3,nslip);
    % S = b(x)m
    % Dyadic product is first-column x second-row vector
    % http://en.wikipedia.org/wiki/Dyadic_product
    for isys=1:nslip,
        s0alpha(:, :, isys) = kron(b0(:,isys), m0(:, isys)');
    end
    % Schmid factor for slip
    % Both +ve and -ve Schmid factors are allowed.
    % Hence the abs()
    PSCHMID(itex) = max(abs(s0alpha(schmid_direction,schmid_direction,:)));
    PSCHMIDALL = squeeze(abs(s0alpha(schmid_direction,schmid_direction,:)));
    % Schmid factor for transformation
    % 3 variants
    ntrans = size(E, 1);
    % get b, m in crystal coordinates
    % [cb0,cm0] = slipsystem_t(nslip, 'T');
    % Transform to sample
    Eg = zeros(3,3,ntrans);
    for i=1:ntrans
        Eg(:, :, i) = tlg*E{i}*tlg';
    end
    s0alpha = Eg;
    % -ve Schmid should loose!
    TVSCHMID(itex) = max(loading_type*s0alpha(schmid_direction,schmid_direction,:));
    TVSCHMIDALL = squeeze(loading_type*s0alpha(schmid_direction,schmid_direction,:))';
    % Schmid factor for plates
    nslip = size(cb0t, 2);
    % Transform to sample
    b0=tlg*cb0t;
    m0=tlg*cm0t;
    s0alpha = zeros(3, 3, nslip);
    for isys=1:nslip,
        s0alpha(:, :, isys) = kron(b0(:,isys), m0(:, isys)');
    end
    % -ve Schmid should loose!
    TPSCHMID(itex) = max(loading_type*s0alpha(schmid_direction,schmid_direction,:));
    TPSCHMIDALL = squeeze(loading_type*s0alpha(schmid_direction,schmid_direction,:))';
    %
end

clear Eg TEXSIZE TEX_POLYand b0 m0 cb0 cm0 euler i isys itex m0;
clear ntrans nslip s0alpha tlg;