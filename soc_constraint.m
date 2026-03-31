function [c,ceq] = soc_constraint(z, z_prev, eps, Nx, Ny, Ngam, N)
        x_prev   = z_prev(Nx    +1:Nx   + N);
        y_prev   = z_prev(Ny    +1:Ny   + N);
        gam_prev = z_prev(Ngam  +1:Ngam + N);
        
        x   = z(Nx    +1:Nx   + N);
        y   = z(Ny    +1:Ny   + N);
        gam = z(Ngam  +1:Ngam + N);
        
        %% SOC of r:
        cx      = norm(x    - x_prev, 2)    - eps;
        cy      = norm(y    - y_prev, 2)    - eps;
        cgam    = norm(gam  -gam_prev, 2)   - eps;
        
        c = [
                cx;
                cy;
                cgam;
             ];
         ceq = [];
         
        
end