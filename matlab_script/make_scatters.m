function res = make_scatters(phantom_img_path, phantom_matrix_path, N)

    img = imread(phantom_img_path);

    img = img';
    [Nl, Ml] = size(img);

    % Define image coordinates


    x_size = 48/1000 ;    %  Size in x-direction [m]
    dx=x_size/Nl;          %  Sampling interval in x direction [m]
    z_size = 50/1000 ;    %  Size in z-direction [m]
    dz=z_size/Ml;          %  Sampling interval in z direction [m]

    y_size = 0/1000;      %  Size in y-direction [m]

    z_start = 5/1000;

    theta = 0;
    % Calculate position data

    x0 = rand(N, 1);
    x = (x0-0.5)* x_size;
    z0 = rand(N, 1);
    z = z0*z_size+z_start;

    y0 = rand(N, 1);
    y = (y0 - 0.5)* y_size;

    %  Find the index for the amplitude value

    xindex = round((x + 0.5*x_size)/dx + 1);
    zindex = round((z - z_start)/dz + 1);
    inside = (0 < xindex)  & (xindex <= Nl) & (0 < zindex)  & (zindex <= Ml);
    index = (xindex + (zindex-1)*Nl).*inside + 1*(1-inside);

    % Amplitudes with different variance must be generated according to the the
    % input map.
    % The amplitude of the liver-kidney image is used to scale the variance

    img_amp = double(img(index) * 1.0)/100.0;

    amp=exp(img_amp);

    amp=amp-min(min(amp));
    amp=1e6*amp/max(max(amp));
    % amp=amp.*randn(N,1).*inside;
    amp = amp.*inside;

    %  Generate the rotated and offset block of sample

    xnew=x*cos(theta)+z*sin(theta);
    znew=z*cos(theta)-x*sin(theta);
    znew=znew-min(min(znew)) + z_start;

    positions=[xnew y znew];

    phantom = us_phantom();
    phantom.sca = positions;
    phantom.amp = amp;
    phantom.N_scatterers = length(positions);


    phantom.write_file(phantom_matrix_path);

    res = 0;

end