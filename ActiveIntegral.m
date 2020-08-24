classdef ActiveIntegral
    
    properties
        Rogowski
        R
        C
        C1
        C2
        R2
        R3
        R4
        R5
        Rf
    end
    
    methods
        function sys = TransferFunction(obj)
             s = tf('s');
             
             R0 = obj.Rogowski.R0;
             L0 = obj.Rogowski.L0;
             C0 = obj.Rogowski.C0;
             Rs = obj.Rogowski.Rs;
             M = obj.Rogowski.M;

             sys = ( ...
                       (M * Rs * obj.R5) / ...
                       ((R0 + Rs) * obj.R4) * s ...
                   ) ...
                   / ...
                   ( ...
                       (L0 * C0 * Rs) / (R0 + Rs) * s * s + ...
                       (L0 + R0 * Rs * C0) / (R0 + Rs) * s + ...
                       1 ...
                   ) ...
                   * ...
                   ( ...
                       (obj.Rf * obj.C1 * s) ...
                       / ...
                       ( ...
                           (obj.Rf * obj.C * s + 1) * ...
                           (obj.R * obj.C1 * s + 1) ...
                       ) ...
                   );
        end
        
        function res = CalcOmegaH(obj)
            R0 = obj.Rogowski.R0;
            Rs = obj.Rogowski.Rs;
            L0 = obj.Rogowski.L0;
            C0 = obj.Rogowski.C0;
            
            res = 0.38 * sqrt((R0 + Rs) / ...
                 (L0 * C0 * Rs));
        end
        
        function res = CalcOmegaL(obj)
             res = 7.07 * 1 / (obj.Rf * obj.C);
        end
        
        function res = CalcSensitivity(obj)
            M = obj.Rogowski.M;
            Rs = obj.Rogowski.Rs;
            R0 = obj.Rogowski.R0;
            res = M * Rs * obj.R5/ ...
                ((R0 + Rs) * obj.R4 * obj.R * obj.C);
        end
    end
end

