function [xee2,Pk2] = Function_MSE_KF_F(F,H,xee2,yy_MSE,Pk2,Q,R)
%MSE_KF_F 此处显示有关此函数的摘要
%   此处显示详细说明
        Pke2 = F * Pk2 * F'+ Q;% ;
        G_MSE = Pke2 * H' * inv(H*Pke2*H' + R);
        xee2 = F * xee2 + G_MSE*(yy_MSE -H*F*xee2);
        Pk2 = Pke2 - G_MSE*H*Pke2;
end

