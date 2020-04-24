classdef JumpPhase < handle
    properties
        t0 = [];
        t1 = [];
        v = [];
        h = [];
        l1 = [];
        l2 = [];
        l3 = [];
        m1 = [];
        m2 = [];
        m3 = [];
        phi = [];
        g = 9.80665;
        vx = []
        vy = 0.0;
        vz = [];
        dq1 = 0.0;
        dq2 = [];
        dq3 = [];
        q1 = [];
        q2 = [];
        q3 = [];
        ddq1 = [];
        ddq2 = [];
        ddq3 = [];
        tau1 =[];
        tau2 = [];
        tau3 = [];
    end
    
    methods
        function self = JumpPhase(t0, t1, v, h, l1, l2, l3, m1, m2, m3, phi)
            self.t0 = t0;
            self.t1 = t1;
            self.v = v;
            self.h = h;
            self.l1 = l1;
            self.l2 = l2;
            self.l3 = l3;
            self.m1 = m1;
            self.m2 = m2;
            self.m3 = m3;
            self.phi = phi;
        end
        
        function vxCoG(self, t)
            beta = self.v*1000/3600/(self.t1-self.t0);
            self.vx = -beta*t;
        end
        
        function vzCoG(self, t)
            % boundary condition
            vz0 = 0.0;
            vz1 = sqrt(2*self.g*self.h);
            dvz0 = 0.0;
            dvz1 = self.g;
            d2vz0 = 0.0;
            d2vz1 = 0.0;
            % spline5
            u = [self.t0, self.t1];
            bound = [vz0, vz1, dvz0, dvz1, d2vz0, d2vz1];
            X = spline5(u, bound);
            self.vz = X(1)*t^5+X(2)*t^4+X(3)*t^3+X(4)*t^2+X(5)*t+X(6);
%             az = 5*X(1)*t^4+4*X(2)*t^3+3*X(3)*t^2+2*X(4)*t+X(5);
%             daz=20*X(1)*t^3+12*X(2)*t^2+6*X(3)+2*X(4);
        end
        
        function angvel(self, t, ~, q2_start, q3_start)
            self.vxCoG(t);
            self.vzCoG(t);
%             V = [self.vx self.vz]';
%             J = [-self.l2*cos(q2_start) -self.l3*cos(q3_start);
%                 -self.l2*sin(q2_start) -self.l3*sin(q3_start)];
            V = [self.vx self.vz]';
            J = [-self.m1*self.l2*cos(q2_start-self.phi)-self.m2*self.l2*cos(q2_start)/2 -(self.m1+self.m2+self.m3/2)*self.l3*cos(q3_start);
                -self.m1*self.l2*sin(q2_start-self.phi)-self.m2*self.l2*sin(q2_start)/2 -(self.m1+self.m2+self.m3/2)*self.l3*sin(q3_start)]...
                /(self.m1+self.m2+self.m3);
            dq = J\V;
            self.dq2 = dq(1);
            self.dq3 = dq(2);
        end
        
        function torque(self, t, q1_start, q2_start, q3_start)
            self.angvel(t, q1_start, q2_start, q3_start);
            self.q1 = q1_start; % self.q1 = int(self.dq1)+q1_start;
            self.q2 = int(self.dq2)+q2_start;
            self.q3 = int(self.dq3)+q3_start;
            self.ddq1 = 0; % self.ddq1 = diff(self.dq1);
            self.ddq2 = diff(self.dq2);
            self.ddq3 = diff(self.dq3);
            self.tau1 = self.m1/2*(self.l1^2*self.ddq1/2+self.l1*self.l2*self.ddq2*cos(self.q1-self.q2+self.phi)+self.l1*self.l2*self.ddq2*sin(self.q1-self.q2+self.phi)+self.l1*self.l3*self.ddq3*cos(self.q1-self.q3)+self.l1*self.l3*self.dq3^2*sin(self.q1-self.q3)+2*self.l2*self.l3*self.dq2*self.dq3*sin(self.q2-self.q3-self.phi)-self.g*self.l1*sin(self.q1));
            self.tau2 = self.m1/2*(2*self.l2^2*self.ddq2+self.l1*self.l2*self.ddq1*cos(self.q1-self.q2+self.phi)-self.l1*self.l2*self.dq1^2*sin(self.q1-self.q2+self.phi)+2*self.l2*self.l3*self.ddq3*cos(self.q2-self.q3-self.phi)+2*self.l2*self.l3*self.dq3^2*sin(self.q2-self.q3-self.phi)-2*self.g*self.l2*sin(self.q2-self.phi))...
                +self.m2/2*(self.l2^2*self.ddq2/2+self.l2*self.l3*self.ddq3*cos(self.q2-self.q3)+self.l2*self.l3*self.dq3^2*sin(self.q2-self.q3)-2*self.g*self.l2*sin(self.q2));
            self.tau3 = self.m1/2*(2*self.l3^2*self.ddq3+self.l1*self.l3*self.ddq1*cos(self.q1-self.q3)-self.l1*self.l3*self.dq1^2*sin(self.q1-self.q3)+2*self.l2*self.l3*self.ddq2*cos(self.q2-self.q3-self.phi)+2*self.l2*self.l3*self.dq2^2*sin(self.q2-self.q3-self.phi)-2*self.g*self.l3*sin(self.q3))...
                +self.m2/2*(2*self.l3^2*self.ddq3+self.l2*self.l3*self.ddq2*cos(self.q2-self.q3)-self.l2*self.l3*self.dq2^2*sin(self.q2-self.q3)-2*self.g*self.l3*sin(self.q3))...
                +self.m3/4*(self.l3^2*self.ddq3-2*self.g*self.l3*sin(self.q3));
        end
    end
end
            
            
            
            
            
        
        
        
        