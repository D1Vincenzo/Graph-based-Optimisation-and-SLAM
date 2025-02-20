classdef LandmarkRangeBearingEdge < g2o.core.BaseBinaryEdge
    % LandmarkRangeBearingEdge summary of LandmarkRangeBearingEdge
    %
    % This class stores an edge which represents the factor for observing
    % the range and bearing of a landmark from the vehicle. Note that the
    % sensor is fixed to the platform.
    %
    % The measurement model is
    %
    %    z_(k+1)=h[x_(k+1)]+w_(k+1)
    %
    % The measurements are r_(k+1) and beta_(k+1) and are given as follows.
    % The sensor is at (lx, ly).
    %
    %    dx = lx - x_(k+1); dy = ly - y_(k+1)
    %
    %    r(k+1) = sqrt(dx^2+dy^2)
    %    beta(k+1) = atan2(dy, dx) - theta_(k+1)
    %
    % The error term
    %    e(x,z) = z(k+1) - h[x(k+1)]
    %
    % However, remember that angle wrapping is required, so you will need
    % to handle this appropriately in compute error.
    %
    % Note this requires estimates from two vertices - x_(k+1) and l_(k+1).
    % Therefore, this inherits from a binary edge. We use the convention
    % that vertex slot 1 contains x_(k+1) and slot 2 contains l_(k+1).
    
    methods(Access = public)
    
        function obj = LandmarkRangeBearingEdge()
            % LandmarkRangeBearingEdge for LandmarkRangeBearingEdge
            %
            % Syntax:
            %   obj = LandmarkRangeBearingEdge(landmark);
            %
            % Description:
            %   Creates an instance of the LandmarkRangeBearingEdge object.
            %   Note we feed in to the constructor the landmark position.
            %   This is to show there is another way to implement this
            %   functionality from the range bearing edge from activity 3.
            %
            % Inputs:
            %   landmark - (2x1 double vector)
            %       The (lx,ly) position of the landmark
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a ObjectGPSMeasurementEdge

            obj = obj@g2o.core.BaseBinaryEdge(2);
        end
        
        function initialEstimate(obj)
            % INITIALESTIMATE Compute the initial estimate of the landmark.
            %
            % Syntax:
            %   obj.initialEstimate();
            %
            % Description:
            %   Compute the initial estimate of the landmark given the
            %   platform pose and observation.
            
            % warning('LandmarkRangeBearingEdge.initialEstimate: implement')
            % 
            % lx = obj.edgeVertices{1}.x(1:2);
            % obj.edgeVertices{2}.setEstimate(lx);

            % Retrieve vehicle pose x_(k+1) = [x, y, psi] (Section A.1)
            % psi & beta (used in MATLAB codes) are the same
            vehiclePose = obj.edgeVertices{1}.estimate;
            % disp(vehiclePose);
            x_k1 = vehiclePose(1);
            y_k1 = vehiclePose(2);
            psi_k1 = vehiclePose(3);
            
            % Retrieve sensor measurements z_(k+1)
            measurement = obj.z;
            r_k1 = measurement(1);  % range
            beta_k1 = measurement(2);  % bearing
            
            % % Print measurement values for debugging
            % z_dash = obj.measurement;
            % disp(z_dash);
            
            % Compute landmark position
            l_x = x_k1 + r_k1 * cos(beta_k1 + psi_k1);
            l_y = y_k1 + r_k1 * sin(beta_k1 + psi_k1);
            
            % Assign the estimated landmark position
            obj.edgeVertices{2}.setEstimate([l_x; l_y]);
            
        end
        
        function computeError(obj)
            % COMPUTEERROR Compute the error for the edge.
            %
            % Syntax:
            %   obj.computeError();
            %
            % Description:
            %   Compute the value of the error, which is the difference
            %   between the predicted and actual range-bearing measurement.

            % warning('LandmarkRangeBearingEdge.computeError: implement')
            % 
            % obj.errorZ = zeros(2, 1);

            x_k1 = obj.edgeVertices{1}.estimate;
            l_k1 = obj.edgeVertices{2}.estimate;
            
            % Compute differences
            dx = l_k1(1) - x_k1(1);
            dy = l_k1(2) - x_k1(2);
            
            % Predicted measurements
            r_pred = sqrt(dx^2 + dy^2);
            beta_pred = atan2(dy, dx) - x_k1(3);
            
            % Error
            obj.errorZ(1) = -obj.z(1) + r_pred;
            obj.errorZ(2) = g2o.stuff.normalize_theta(-obj.z(2) + beta_pred);
            
        end
        
        function linearizeOplus(obj)
            % linearizeOplus Compute the Jacobian of the error in the edge.
            %
            % Syntax:
            %   obj.linearizeOplus();
            %
            % Description:
            %   Compute the Jacobian of the error function with respect to
            %   the vertex.
            %

            % warning('LandmarkRangeBearingEdge.linearizeOplus: implement')
            % 
            % obj.J{1} = eye(2, 3);
            % 
            % obj.J{2} = eye(2);

            x_k1 = obj.edgeVertices{1}.estimate;
            l_k1 = obj.edgeVertices{2}.estimate;
            
            % Compute differences
            dx = l_k1(1) - x_k1(1);
            dy = l_k1(2) - x_k1(2);
            
            % Predicted range
            r_pred = sqrt(dx^2 + dy^2);
            r_sq = max(r_pred^2, 1e-6);
            
            % Compute Jacobians
            J1 = zeros(2, 3); % with respect to x_k+1
            J2 = zeros(2, 2); % with respect to l_k+1
            
            J1(1, 1) = -dx / r_pred;
            J1(1, 2) = -dy / r_pred;
            J1(1, 3) = 0;
            
            J1(2, 1) = dy / r_sq;
            J1(2, 2) = -dx / r_sq;
            J1(2, 3) = -1;
            
            J2(1, 1) = dx / r_pred;
            J2(1, 2) = dy / r_pred;
            
            J2(2, 1) = -dy / r_sq;
            J2(2, 2) = dx / r_sq;
            
            % Assign Jacobians
            obj.J{1} = J1;
            obj.J{2} = J2;

        end        
    end
end