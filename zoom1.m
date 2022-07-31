function alphas = zoom1(f0,f1,x0,d,alphal,alphah)
% function alphas = zoom(f,x0,d,alphal,alphah)
% Algorithm 3.6 on page 61 in Nocedal and Wright

c1 = 0.25;
c2 = 0.75;
fx0 = f0(x0(1),x0(2));
gx0 = f1(x0(1),x0(2))'*d;

while (1~=2),
   alphax = 1/2*(alphal+alphah);
   xx = x0 + d.*alphax;
   fxx = f0(xx(1),xx(2));
   gxx = f1(xx(1),xx(2))'*d;
   xl = x0 + d.*alphal;
   fxl = f0(xl(1),xl(2));
   gxl = f1(xl(1),xl(2));
   if ((fxx > fx0 + c1*alphax*gx0) || (fxx >= fxl)),
      alphah = alphax;
   else
      if abs(gxx) < -c2*gx0,
        alphas = alphax;
        return;
      end
      if gxx*(alphah-alphal) >= 0,
        alphah = alphal;
      end
      alphal = alphax;
   end
end 