classdef scan < uff
    %SCAN   UFF class to define a scan 
    %   SCAN contains the position of a collection of pixels. It is a
    %   superclass for more easy-to-handle classes such as UFF.LINEAR_SCAN
    %   or UFF.SECTOR_SCAN
    %
    %   Compulsory properties:
    %         x                  % Vector containing the x coordinates of each pixel in [m]
    %         y                  % Vector containing the y coordinates of each pixel in [m]
    %         z                  % Vector containing the z coordinates of each pixel in [m]
    %
    %   Example:
    %         sca = uff.scan();
    %         x_axis=linspace(-20e-3,20e-3,256);
    %         z_axis=linspace(0e-3,40e-3,256);
    %         [X Z]=meshgrid(x_axis,z_axis);
    %         sca.x=X(:);
    %         sca.y=zeros(size(X(:)));
    %         sca.z=Z(:);
    %
    %   See also UFF.LINEAR_SCAN, UFF.SECTOR_SCAN


    properties  (Access = public)
        x                  % Vector containing the x coordinate of each pixel in the matrix
        y                  % Vector containing the x coordinate of each pixel in the matrix
        z                  % Vector containing the z coordinate of each pixel in the matrix
    end
    
    properties  (Dependent)
        N_pixels           % total number of pixels in the matrix
        xyz                % location of the source [m m m] if the source is not at infinity    
    end
    
    %% constructor -> uff constructor
    methods (Access = public)
        function h=scan(varargin)
            h = h@uff(varargin{:});
        end
    end
    
    %% plot methods
    methods
        function figure_handle=plot(h,figure_handle_in,title_in)
            % plotting scan
            if (nargin>1) && ~isempty(figure_handle_in)
                figure_handle=figure(figure_handle_in); hold on;
            else
                figure_handle=figure();
                title('Probe');
            end
            
            plot3(h.x*1e3,h.y*1e3,h.z*1e3,'k.');
            xlabel('x[mm]'); ylabel('y[mm]'); zlabel('z[mm]');
            set(gca,'ZDir','Reverse');
            set(gca,'fontsize',14);
            
            if nargin>2
                title(title_in);
            end
        end
    end
    
    %% Set methods
    methods
        function h=set.x(h,in_x)
            assert(size(in_x,2)==1, 'The x vector must be a column vector.')
            h.x=in_x;
        end
        function h=set.y(h,in_y)
            assert(size(in_y,2)==1, 'The y vector must be a column vector.')
            h.y=in_y;
        end
        function h=set.z(h,in_z)
            assert(size(in_z,2)==1, 'The z vector must be a column vector.')
            h.z=in_z;
        end
        function h=set.xyz(h,in_xyz)
             assert(size(in_xyz,2)==3, 'The xyz must be an array [x y z] - [m m m]');
             h.x=in_xyz(:,1);
             h.y=in_xyz(:,2);
             h.z=in_xyz(:,3);
        end
    end
    
    %% Get methods
    methods
        function value=get.N_pixels(h)
            value=min([numel(h.x) numel(h.y) numel(h.z)]);
        end
        function value=get.xyz(h)
             value=[h.x h.y h.z];
        end
    end
end

