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

            % Compute the posterior assming no noise
            x_k = obj.edgeVertices{1}.x; % current stage (x_k=[x_k, y_k, psi_k])
            psi_k = x_k(3);
            u_k = obj.measurement; % The control input
            
            if isfield(obj.edgeVertices{1}, 'time') && isfield(obj.edgeVertices{2}, 'time')
                deltaT = abs(obj.edgeVertices{2}.time - obj.edgeVertices{1}.time);
            else
                deltaT = obj.dT;
            end


            M = deltaT*[cos(psi_k), -sin(psi_k), 0;
                        sin(psi_k),  cos(psi_k), 0;
                        0, 0, 1];

            x_k1 = x_k + M * u_k;  % assming no noise

            % Assign computed state to the next vertex
            obj.edgeVertices{2}.x = x_k1; 

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
            
            x_k = obj.edgeVertices{1}.x; % current stage (x_k=[x_k, y_k, psi_k])
            psi_k = x_k(3);
            x_k1 = obj.edgeVertices{2}.x;  % Next state from graph
            u_k = obj.z;

            if isfield(obj.edgeVertices{1}, 'time') && isfield(obj.edgeVertices{2}, 'time')
                deltaT = abs(obj.edgeVertices{2}.time - obj.edgeVertices{1}.time);
            else
                deltaT = obj.dT;
            end

            M = deltaT*[cos(psi_k), -sin(psi_k), 0;
                        sin(psi_k),  cos(psi_k), 0;
                        0, 0, 1];

            det_M = det(M); % 計算行列式
            disp(['Determinant of M: ', num2str(det_M)]);

            % obj.errorZ = pinv(M) \ (x_k1 - x_k);
            % obj.errorZ = (M \ (x_k1 - x_k)) - u_k;

            x_k = x_k + M * u_k;
            obj.errorZ = M \ (x_k1 - x_k);

            % x_k1 = obj.edgeVertices{2}.x;  % Next state from graph
            % 
            % predicted_x_k1 = obj.edgeVertices{2}.x;  % predict next stage with initialEstimate()
            % obj.errorZ = x_k1 - predicted_x_k1; 


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

            % Since computeError() now applies inv(M), J_x should just be -I
            J_x = -eye(3);
            
            % Jacobian w.r.t. x_k1 is identity since x_k1 directly contributes to itself
            J_x1 = eye(3);
            
            % Store Jacobians
            obj.J{1} = J_x;  % Influence of previous state
            obj.J{2} = J_x1; % Influence of next state
        end
    end    
end