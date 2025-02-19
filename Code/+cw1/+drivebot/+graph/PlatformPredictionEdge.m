classdef PlatformPredictionEdge < g2o.core.BaseBinaryEdge
    % PlatformPredictionEdge summary of PlatformPredictionEdge
    %
    % This class stores the factor representing the process model which
    % transforms the state from timestep k to k+1
    %
    % The process model is as follows.
    %
    % Define the rotation vector
    %
    %   M = dT * [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0;0 0 1];
    %
    % The new state is predicted from 
    %
    %   x_(k+1) = x_(k) + M * [vx;vy;theta]
    %
    % Note in this case the measurement is actually the mean of the process
    % noise. It has a value of 0. The error vector is given by
    %
    % e(x,z) = inv(M) * (x_(k+1) - x_(k))
    %
    % Note this requires estimates from two vertices - x_(k) and x_(k+1).
    % Therefore, this inherits from a binary edge. We use the convention
    % that vertex slot 1 contains x_(k) and slot 2 contains x_(k+1).
    
    properties(Access = protected)
        % The length of the time step
        dT;
    end
    
    methods(Access = public)
        function obj = PlatformPredictionEdge(dT)
            % PlatformPredictionEdge for PlatformPredictionEdge
            %
            % Syntax:
            %   obj = PlatformPredictionEdge(dT);
            %
            % Description:
            %   Creates an instance of the PlatformPredictionEdge object.
            %   This predicts the state from one timestep to the next. The
            %   length of the prediction interval is dT.
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a PlatformPredictionEdge

            assert(dT >= 0);
            obj = obj@g2o.core.BaseBinaryEdge(3);            
            obj.dT = dT;
        end
       
        function initialEstimate(obj)
            % INITIALESTIMATE Compute the initial estimate of a platform.
            %
            % Syntax:
            %   obj.initialEstimate();
            %
            % Description:
            %   Compute the initial estimate of the platform x_(k+1) given
            %   an estimate of the platform at time x_(k) and the control
            %   input u_(k+1)

            %%%%%
            %%%%%
            % Retrieve the current state estimate from vertex 1 (x_k)
            xk = obj.edgeVertices{1}.x;
            
            % Retrieve the control input from the edge measurement
            if isempty(obj.measurement) || length(obj.measurement) ~= 3
                error('Control input u is missing or has an incorrect size. Expected [vx; vy; dtheta].');
            end
            u = obj.measurement;
            
            % Extract the current heading angle from x_k (assumed to be the third element)
            theta = xk(3);
            
            % Normalize theta to the range [-pi, pi]
            theta = mod(theta + pi, 2*pi) - pi;
            
            % Ensure dT is a valid time step
            dT = obj.dT;
            if isempty(dT) || abs(dT) < 1e-6
                warning('dT is too small or undefined, setting it to 1e-6 to avoid numerical issues.');
                dT = 1e-6;  % Small regularization term
            end
            
            % Build the transformation matrix M
            M = dT * [ cos(theta), -sin(theta), 0;
                       sin(theta),  cos(theta), 0;
                       0,           0,          1 ];
            
            % Compute the predicted state at time k+1
            xk1_pred = xk + M * u;
            
            % Ensure valid assignment to vertex 2
            if ~isfield(obj.edgeVertices{2}, 'x')
                warning('Vertex 2 does not have a state variable "x". Creating it.');
                obj.edgeVertices{2}.x = zeros(size(xk));
            end
            
            % Set the initial estimate for vertex 2
            obj.edgeVertices{2}.x = xk1_pred;
            %%%%%
            %%%%%

            % warning('PlatformPredictionEdge.initialEstimate: implement')
            % 
            % % Compute the posterior assming no noise
            % obj.edgeVertices{2}.x = zeros(3, 1);
        end
        
        function computeError(obj)
            % COMPUTEERROR Compute the error for the edge.
            %
            % Syntax:
            %   obj.computeError();
            %
            % Description:
            %   Compute the value of the error, which is the difference
            %   between the measurement and the parameter state in the
            %   vertex. Note the error enters in a nonlinear manner, so the
            %   equation has to be rearranged to make the error the subject
            %   of the formulat

            %%%%%
            %%%%%
            % Retrieve the current estimates from the two vertices.
            xk   = obj.vertex(1).estimate();  % x_(k)
            xk1  = obj.vertex(2).estimate();  % x_(k+1)
            
            % Extract the heading angle from x_(k).
            theta = xk(3);
            
            % Normalize theta to the range [-pi, pi]
            theta = mod(theta + pi, 2*pi) - pi;
            
            % Avoid division by zero by ensuring dT is nonzero
            dT = obj.dT;
            if abs(dT) < 1e-6
                warning('dT is too small, adding regularization to avoid singularity.');
                dT = 1e-6;  % Small regularization term
            end
            
            % Construct the matrix M based on the heading and the time step dT.
            M = dT * [ cos(theta), -sin(theta), 0; 
                       sin(theta),  cos(theta), 0;
                       0,           0,          1 ];
            
            % Compute the error: e = inv(M) * ( x_(k+1) - x_(k) ).
            % Using left-division to avoid explicit inversion.
            obj.errorZ = M \ (xk1 - xk);
            %%%%%
            %%%%%

            % warning('PlatformPredictionEdge.computeError: implement')

            % obj.errorZ = 0;
        end
        
        % Compute the Jacobians
        function linearizeOplus(obj)
            % LINEARIZEOPLUS Compute the Jacobians for the edge.
            %
            % Syntax:
            %   obj.computeError();
            %
            % Description:
            %   Compute the Jacobians for the edge. Since we have two
            %   vertices which contribute to the edge, the Jacobians with
            %   respect to both of them must be computed.
            %

            %%%%%
            %%%%%
            % Retrieve current estimates from the two vertices:
            if isempty(obj.edgeVertices{1}.x) || isempty(obj.edgeVertices{2}.x)
                error('State estimates for vertices are missing.');
            end
            
            xk  = obj.edgeVertices{1}.x;  % State at time k: [x; y; theta]
            xk1 = obj.edgeVertices{2}.x;  % State at time k+1
        
            % Extract the heading angle theta from xk:
            theta = xk(3);
            
            % Normalize theta to range [-pi, pi]
            theta = mod(theta + pi, 2*pi) - pi;
            
            % Ensure dT is valid
            dT = obj.dT;
            if isempty(dT) || abs(dT) < 1e-6
                warning('dT is too small or undefined. Setting dT to 1e-6 for stability.');
                dT = 1e-6;
            end
        
            % Build the rotation matrix R (for a 2D rotation) extended to 3x3:
            R = [ cos(theta), -sin(theta), 0;
                  sin(theta),  cos(theta), 0;
                  0,           0,          1 ];
              
            % Construct the process model matrix M = dT * R.
            M = dT * R;
            
            % Compute its inverse: since M = dT * R, we have inv(M) = R' / dT.
            A = (R') / dT;  % A = inv(M)
            
            % Compute the state difference:
            v = xk1 - xk;
            
            % --- Jacobian with respect to vertex 2 ---
            % The error is e = A * (xk1 - xk), so the derivative with respect to xk1 is:
            J2 = A;
            
            % --- Jacobian with respect to vertex 1 ---
            % There are two contributions:
            % 1. From the subtraction: d/dxk ( - (xk1 - xk) ) = -I.
            % 2. From the dependence of A on theta in xk.
            
            % Compute derivative of R' with respect to theta:
            dRtrans_dtheta = [ -sin(theta),  cos(theta), 0;
                              -cos(theta), -sin(theta), 0;
                               0,           0,          0 ];
            
            % Compute the derivative of A = R'/dT with respect to theta:
            dA_dtheta = dRtrans_dtheta / dT;
            
            % The derivative with respect to theta:
            %    d(e)/d(theta) = dA/dtheta * (xk1 - xk)  +  A * (d/dtheta (-xk))
            %                   = dA_dtheta * v - A * [0; 0; 1]
            
            % Construct Jacobian for vertex 1:
            J1 = zeros(3,3);
            J1(:,1:2) = -A(:,1:2);  % First two columns: simple subtraction effect
            J1(:,3) = dA_dtheta * v - A * [0; 0; 1];  % Third column: theta dependency
            
            % Assign computed Jacobians:
            obj.J{1} = J1;
            obj.J{2} = J2;
            %%%%%
            %%%%%

            % warning('PlatformPredictionEdge.linearizeOplus: implement')
            % 
            % obj.J{1} = -eye(3);
            % 
            % obj.J{2} = eye(3);
        end
    end    
end