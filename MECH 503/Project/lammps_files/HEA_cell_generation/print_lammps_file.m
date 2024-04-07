function [output] = print_lammps_file(atoms, Pars, MyRandomGeneration, generations, filename)

output = 0;

xcoord=atoms(:,1);
ycoord=atoms(:,2);
zcoord=atoms(:,3);

a0 = Pars.latconst;
Nx = Pars.Nx;
Ny = Pars.Ny;
Nz = Pars.Nz;

xlo = 0;     %min(xcoord);
xhi = Nx*a0; %max(xcoord);

ylo = 0;     %min(xcoord);
yhi = Ny*a0; %max(xcoord);

zlo = 0;     %min(xcoord);
zhi = Nz*a0; %max(xcoord);

Na = Pars.Na; 
Ne = Pars.Ne; 

k= MyRandomGeneration;

fileID = fopen(filename,'w');

% Print header for lammps
fprintf(fileID,'# File generated with cnn_vac_diff.m by M. Ponga\n\n');
fprintf(fileID,'\t %d  atoms\n',Na);
fprintf(fileID,'\t \t %d  atom types\n\n',Ne);
fprintf(fileID,'\t %f \t %f  xlo xhi\n', xlo, xhi);
fprintf(fileID,'\t %f \t %f  ylo yhi\n', ylo, yhi);
fprintf(fileID,'\t %f \t %f  zlo zhi\n\n', zlo, zhi);
fprintf(fileID,' Masses\n\n');
fprintf(fileID,'\t 1   55.84500000    # Fe\n');
fprintf(fileID,'\t 2   58.69340000    # Ni\n');
fprintf(fileID,'\t 3   51.99610000    # Cr\n');
fprintf(fileID,'\t 4   58.93319500    # Co\n\n');
fprintf(fileID,' Atoms # atomic\n\n');

%
for i=1:Na
    fprintf(fileID,'\t %d \t %d \t %f \t %f \t %f\n', i, generations(k,i), xcoord(i)*a0, ycoord(i)*a0, zcoord(i)*a0);
end

fclose(fileID);

output = 1;

end

