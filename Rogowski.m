classdef Rogowski < handle
    properties
        d % 内径
        D % 外径
        h % 高度
        a % 布线宽度
        N % 匝数
        Cu % 铜箔厚度
        l % 匝间距离
        Resistivity % 电阻率
        
        R0
        L0
        C0
        Rs
        M
    end
    
    methods
       function obj = CalcElecData(obj)
           u0 = 1.2566370614e-6;
           vd0 = 8.854187817e-12;
           
           c = (obj.D - obj.d) / 2;
           
           obj.M = (obj.N * u0 * obj.h)...
               /...
               (2 * pi) * ...
               log(obj.D/obj.d);
           
           obj.L0 = (obj.N * obj.N * u0 * obj.h)...
               /...
               (2 * pi) * ...
               log(obj.D/obj.d);
           
           obj.R0 = (2 * obj.Resistivity * (c + obj.h) * obj.N)...
               /...
               (obj.a * obj.Cu);

           S = 2 * (c + obj.h) * obj.Cu;
           
           obj.C0 = vd0 * 4.5 * S / (obj.l * (obj.N - 1));
       end
       
       function rs = CalcBestRs(obj)
           rs = obj.L0 / sqrt(2 * obj.L0 * obj.C0 - ...
                obj.R0 ^ 2 * obj.C0 ^ 2);
       end
    end
end
