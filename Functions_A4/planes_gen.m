function Planes = planes_gen()

strucfields = struct('plane', [], ... % Plane number
                     'Tf',    [], ... % T_final
                     'umax',  [], ... % u_max
                     'A',     [], ... % A matrix
                     'B',     [], ... % B matrix
                     'dim',   [], ... % Dimensions
                     'x0',    [], ... % initial cond.
                     'T',     [], ... % State pred.
                     'S',     [], ... % Input pred.
                     'H',     [], ... % Cost fun H
                     'h',     [],  ... % Cost fun h
                     'A_eq',  [], ... % x_f constraints
                     'b_eq',  [], ... % x_f constraints
                     'A_u',   [], ... % Control input constraints
                     'b_u',   []  ... % Control input constraints
                    );

P1 = strucfields;
P2 = strucfields;
P3 = strucfields;
P4 = strucfields;
Planes = [P1, P2, P3, P4];

end