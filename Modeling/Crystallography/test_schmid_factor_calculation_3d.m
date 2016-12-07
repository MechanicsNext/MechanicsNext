% Test Schmid factor calculation

% NiTi
%texture_ip          = [9.329181291088404215e-01	1.131258077512668075e-01	-3.323692948023018179e-01	7.998104696304578209e-02]; % 51NiTiSC0D4 cubic orientation
texture_ip          = [1	0	0	0]; % 51NiTiSC0D4 cubic orientation
schmid_direction_ip = 2; % Y loading
loading_type_ip     = 1; % tension
% Transformation matrices (U) for cubic -> monoclinic in NiTi
alpha = 1.0243;
gamma = 0.9563;
delta = 0.058;
epsi = -0.0427;
U_ip = { [gamma  epsi   epsi; % 1
    epsi   alpha  delta;
    epsi   delta  alpha];
    [gamma -epsi  -epsi;      % 2
    -epsi   alpha  delta;
    -epsi   delta  alpha];
    [gamma -epsi   epsi;      % 3
    -epsi   alpha -delta;
    epsi  -delta  alpha];
    [gamma  epsi  -epsi;      % 4
    epsi   alpha -delta;
    -epsi  -delta  alpha];
    [alpha  epsi   delta;     % 5
    epsi   gamma  epsi;
    delta  epsi  alpha];
    [alpha -epsi   delta;     % 6
    -epsi   gamma -epsi
    delta -epsi  alpha];
    [alpha -epsi  -delta;     % 7
    -epsi   gamma  epsi;
    -delta  epsi  alpha];
    [alpha  epsi  -delta;     % 8
    epsi   gamma -epsi;
    -delta -epsi  alpha];
    [alpha  delta  epsi;      % 9
    delta  alpha  epsi;
    epsi   epsi   gamma];
    [alpha  delta -epsi;      % 10
    delta  alpha -epsi;
    -epsi  -epsi   gamma];
    [alpha -delta  epsi;      % 11
    -delta  alpha -epsi;
    epsi  -epsi   gamma];
    [alpha -delta -epsi;      % 12
    -delta  alpha  epsi;
    -epsi   epsi   gamma];
    };
%
[TPSCHMID, TPSCHMIDALL, TVSCHMID, TVSCHMIDALL, PSCHMID, PSCHMIDALL] = schmid_factor_calculation_3d(texture_ip, schmid_direction_ip, loading_type_ip, U_ip);