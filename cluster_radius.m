function [r_angle_deg,r_circle] = cluster_radius(points, clusters, centers, varargin)
%CLUSTER_ANGULAR_RADIUS Per-cluster angular radius for broken/irregular PFs.
%
% points  : (N x M) objective vectors
% idx     : (N x 1) integer labels in {1..K}
% centers : (K x M) cluster centers in the same label order
%
% Optional name-value pairs:
%   'SingleMarginDeg' (default 1.0) -> subtract this from single-member radii
%   'Eps'             (default 1e-12) -> numerical epsilon for norms/cosines
%
% r_angle_deg : (K x 1) angular radius (degrees) per cluster
%
% Notes
% - Angles are measured at the origin: angle(p, c_k).
% - For single-member clusters, midpoints m_{k,j} = (c_k + c_j)/2 define
%   “bisector-like” directions separating neighboring subregions; we take
%   the *smallest* angle between c_k and any m_{k,j} as the boundary.
% - A small margin (SingleMarginDeg) is subtracted to avoid overlap.

p = inputParser;
p.addParameter('SingleMarginDeg', 1.0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('Eps', 1e-12, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.parse(varargin{:});
marginDeg = p.Results.SingleMarginDeg;
epsval    = 1e-12;

K = size(centers,1);
r_angle_deg = zeros(K,1);

% Precompute norms once
point_norms  = sqrt(sum(points.^2, 2));
center_norms = sqrt(sum(centers.^2, 2));

r_angle_deg  = zeros(K,1);
r_circle     = zeros(K,1);

if size(centers,2)==2
            if ~exist('epsval', 'var')
                epsval = 1e-12; 
            end
            if length(clusters)>1
                K = size(centers, 1);
                r_circle     = zeros(K, 1);
                r_angle_deg  = zeros(K, 1);

                % Step 1: Sequential midpoints between consecutive centers
                midpoints = (centers(1:end-1,:) + centers(2:end,:)) / 2;

                % Step 2: Norms
                center_norms   = sqrt(sum(centers.^2, 2));
                midpoint_norms = sqrt(sum(midpoints.^2, 2));

                % Step 3: For each center, choose the best of its adjacent midpoints
                for k = 1:K
                    ck = centers(k,:);
                    nc = max(center_norms(k), epsval);

                    adj_mp = [];
                    adj_idx = [];

                    if k > 1
                        adj_mp = [adj_mp; midpoints(k-1,:)];
                        adj_idx = [adj_idx; k-1];
                    end
                    if k < K
                        adj_mp = [adj_mp; midpoints(k,:)];
                        adj_idx = [adj_idx; k];
                    end

                    best_angle = inf;
                    best_radius = 0;

                    for i = 1:size(adj_mp,1)
                        mk = adj_mp(i,:);
                        nm = max(norm(mk), epsval);

                        cos_theta = dot(ck, mk) / (nc * nm);
                        cos_theta = min(1, max(-1, cos_theta));  % safety clamp
                        angle_deg = acosd(cos_theta);
                        radius    = norm(ck - mk);

                        if angle_deg < best_angle
                            best_angle = angle_deg;
                            best_radius = radius;
                        end
                    end

                    r_angle_deg(k) = best_angle;
                    r_circle(k)    = best_radius;
                end
                for k = 1:K
                    members = [clusters(k).sol];
                    ck  = centers(k,:);
                    nc  = max(center_norms(k), epsval);

                    if isempty(members) || nc < epsval
                        r_angle_deg(k) = 0;
                        r_circle1(k)    = 0;
                        continue;
                    end

                    if numel(members) >= 1

                        pm   = points(members,:);
                        pn   = max(point_norms(members), 1e-12);
                        cos_t = (pm*ck.') ./ (pn*nc);
                        cos_t = max(-1, min(1, cos_t));
                        angles = acosd(cos_t);
                        r_angle_deg1(k) = max(angles) + 1e-6;  % pad a hair

                        % --- radial cap from members (or use prctile to resist outliers) ---
                        % angular radius from your code
                        th = deg2rad(r_angle_deg1(k));

                        % radial band of the points you want guaranteed inside
                        pn_k   = sqrt(sum(pm.^2,2));
                        rmin_k = min(pn_k);
                        rmax_k = max(pn_k);     % or prctile(pn_k,95) to ignore outliers

                        % distances needed at the two endpoints (no origin)
                        d_min = sqrt(max(0, nc^2 + rmin_k^2 - 2*nc*rmin_k*cos(th)));
                        d_max = sqrt(max(0, nc^2 + rmax_k^2 - 2*nc*rmax_k*cos(th)));

                        % final circle radius (smallest that covers the truncated wedge)
                        r_circle1(k) = max(d_min, d_max);
                    end
                    
                end
                a=find(r_angle_deg1>r_angle_deg');
                if length(a)>=1
                     r_angle_deg(a)=r_angle_deg1(a);
                end
                a=find(r_circle1>r_circle');
                if length(a)>=1
                     r_circle(a)=r_circle1(a);
                end
                

        else
            point_norms  = sqrt(sum(points.^2, 2));
            center_norms = sqrt(sum(centers.^2, 2));
            members = [clusters(1).sol];
            ck   = centers;
            nc   = center_norms;
            pm  = points(members,:);
            pn  = max(point_norms(members), epsval);            % avoid /0
            dots = pm * ck';
            cos_t = dots ./ (pn * nc);
            cos_t = min(1, max(-1, cos_t));                     % numeric safety
            angles = acosd(cos_t);
            r_angle_deg = max(angles);   
             % angular radius from your code
            th = deg2rad(r_angle_deg);

            % radial band of the points you want guaranteed inside
            pn_k   = sqrt(sum(pm.^2,2));
            rmin_k = min(pn_k);
            rmax_k = max(pn_k);     % or prctile(pn_k,95) to ignore outliers

            % distances needed at the two endpoints (no origin)
            d_min = sqrt(max(0, nc^2 + rmin_k^2 - 2*nc*rmin_k*cos(th)));
            d_max = sqrt(max(0, nc^2 + rmax_k^2 - 2*nc*rmax_k*cos(th)));

            % final circle radius (smallest that covers the truncated wedge)
            r_circle = max(d_min, d_max);
        end
else
    if length(clusters)>1

        for k = 1:K
            members = [clusters(k).sol];
            ck  = centers(k,:);
            nc  = max(center_norms(k), epsval);

            if isempty(members) || nc < epsval
                r_angle_deg(k) = 0;
                r_circle(k)    = 0;
                continue;
            end

            if numel(members) >= 10

                pm   = points(members,:);
                pn   = max(point_norms(members), 1e-12);
                cos_t = (pm*ck.') ./ (pn*nc);
                cos_t = max(-1, min(1, cos_t));
                angles = acosd(cos_t);
                r_angle_deg(k) = max(angles) + 1e-6;  % pad a hair

                % --- radial cap from members (or use prctile to resist outliers) ---
                % angular radius from your code
                th = deg2rad(r_angle_deg(k));

                % radial band of the points you want guaranteed inside
                pn_k   = sqrt(sum(pm.^2,2));
                rmin_k = min(pn_k);
                rmax_k = max(pn_k);     % or prctile(pn_k,95) to ignore outliers

                % distances needed at the two endpoints (no origin)
                d_min = sqrt(max(0, nc^2 + rmin_k^2 - 2*nc*rmin_k*cos(th)));
                d_max = sqrt(max(0, nc^2 + rmax_k^2 - 2*nc*rmax_k*cos(th)));

                % final circle radius (smallest that covers the truncated wedge)
                r_circle(k) = max(d_min, d_max);
            else
                % ------ Single-member cluster: use midpoint directions ------
                if K == 1
                    % Only one cluster exists; define a conservative radius
                    r_angle_deg(k) = 0;  % nothing else to separate from
                    continue;
                end
                pm=centers(setdiff(1:K,k),:);
                pn   = max(center_norms(setdiff(1:K,k)), 1e-12);
                cos_t = (pm*ck.') ./ (pn*nc);
                cos_t = max(-1, min(1, cos_t));
                angles = acosd(cos_t);
                r_angle_deg(k) = min(angles)/2 + 1e-6;  % pad a hair
                th = deg2rad(r_angle_deg(k));

                % radial band of the points you want guaranteed inside
                pn_k   = sqrt(sum(pm.^2,2));
                rmin_k = min(pn_k);
                rmax_k = max(pn_k);     % or prctile(pn_k,95) to ignore outliers

                % distances needed at the two endpoints (no origin)
                d_min = sqrt(max(0, nc^2 + rmin_k^2 - 2*nc*rmin_k*cos(th)));
                d_max = sqrt(max(0, nc^2 + rmax_k^2 - 2*nc*rmax_k*cos(th)));

                % final circle radius (smallest that covers the truncated wedge)
                r_circle(k) = max(d_min, d_max);
            

            end
        end
    else
        point_norms  = sqrt(sum(points.^2, 2));
        center_norms = sqrt(sum(centers.^2, 2));
        members = [clusters(1).sol];
        ck   = centers;
        nc   = center_norms;
        pm  = points(members,:);
        pn  = max(point_norms(members), epsval);            % avoid /0
        dots = pm * ck';
        cos_t = dots ./ (pn * nc);
        cos_t = min(1, max(-1, cos_t));                     % numeric safety
        angles = acosd(cos_t);
        r_angle_deg = max(angles);  
        nc  = norm(centers);
        th  = deg2rad(r_angle_deg);
        % choose the radial cap you want to guarantee
        pn_all = sqrt(sum(points(members,:).^2,2));
        Rk  = max(pn_all);              % or: max norm of cluster k’s members
        r_circle = max(nc, sqrt(nc^2 + Rk^2 - 2*nc*Rk*cos(th)));

    end
end
end

