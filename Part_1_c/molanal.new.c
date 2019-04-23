/********************************************************************

   A simple molecular analyzer for DFTB simulations.

   usage: molanal.new <gen file> <chr file>

   The chr (charge) file is optional.
   The gen (coordinate) file is mandatory.

   Non-orthorhombic simulation cells are supported.

   Output:  Molecules found are printed out with total charges,
         if supplied.
         An xyz file suitable for jmol is created in the
   file molanal.xyz.

   Bond distances are now read from the file bonds.dat.

   Larry Fried, 5/30/2005.

   A bond lifetime criterion is now implemented.  Fixed a bug in
   dipole calculations when wrap_com was not called.

   Larry Fried, 11/10/06

   WARNING: THIS FILE HAS BEEN CUSTOM MODIFIED - Qian Yang, 1/29/2014
   It is designed to work with the ATOM dump file from LAMMPS, and
   furthermore assumes elements are C and H.

   - Reactions are now detected and printed
   - Feature vectors are outputted
   Cameron Jones, 7/3/18

   0.012 picoseconds / frame



 *********************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#define BUFSIZE 1024  /** Size of buffer for IO. **/
#define MAXATOM 1100 /*1100*/  /** Maximum number of atoms. Default: 256 **/
#define MAXWRAP 2   /** Maximum number of lattice vectors to search over.
	                          Use a larger number for a highly non-orthogonal cell. **/
#define MAXELE 4    /** Maximum number of elements **/
#define MAXHIST 21  /** Maximum number of history steps for bond determination. **/
#define NOTSET 9999   /** Number to identify whether history variables have been set. **/

// Define character codes for colored terminal output
#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define MAG   "\x1B[35m"
#define CYN   "\x1B[36m"
#define WHT   "\x1B[37m"
#define RESET "\x1B[0m"

int not_in_molecule(int isearch, int mol_list[MAXATOM][MAXATOM],
                    int nmolecule);
void find_molecules(double *x, double *y, double *z,
                    double a[3], double b[3], double c[3],
                    int mol_list[MAXATOM][MAXATOM],
                    int natom, int *nmolecule,
                    int *type, char *element,
                    int bond_list[MAXATOM][MAXATOM]);
int is_bonded(int j, int k, double *x, double *y, double *z,
              double a[3], double b[3], double c[3],
              int *type, char *element, int mol_list[MAXATOM][MAXATOM], int nmolecule);
void find_molecule(int j,
                   double *x, double *y, double *z,
                   double a[3],  double b[3], double c[3],
                   int mol_list[MAXATOM][MAXATOM],
                   int natom,
                   int nmolecule,
                   int *type,
                   char *element,
                   int bond_list[MAXATOM][MAXATOM]);
void wrap_atoms(double *x, double *y, double *z,
                double a[3], double b[3], double c[3],
                int bond_list[MAXATOM][MAXATOM],
                int natom);
double wrap_atom(int j, int k, double *x, double *y, double *z,
                 double a[3], double b[3], double c[3], int maxwrap,
                 int do_wrap);
void wrap_molecules(double *x, double *y, double *z,
                    double a[3], double b[3], double c[3],
                    int mol_list[MAXATOM][MAXATOM], int nmolecule,
                    int bond_list[MAXATOM][MAXATOM], int natom );
void wrap_com(double *x, double *y, double *z,
              double a[3], double b[3], double c[3],
              int mol_list[MAXATOM][MAXATOM], int nmolecule, int natom,
              double rcomx[MAXATOM], double rcomy[MAXATOM],
              double rcomz[MAXATOM], int type[MAXATOM],
              char element[MAXELE]);
void center_all(double *x, double *y, double *z,
                double a[3], double b[3], double c[3],
                int mol_list[MAXATOM][MAXATOM], int nmolecule);
double r2bond(char ea, char eb);
int nint(double num);

//print_molecule is not used anymore
//void print_molecule(int mol_list[MAXATOM][MAXATOM],
//                    int bond_list[MAXATOM][MAXATOM],
//                    int *type, char *element, int mol,
//                    int read_charge,
//                    double q[MAXATOM],
//                    double x[MAXATOM], double y[MAXATOM], double z[MAXATOM],
//                    double rcomx[MAXATOM], double rcomy[MAXATOM],
//                    double rcomz[MAXATOM], int printed_atom[MAXATOM],
//                    int dump_structures, int nmolecule, int real_frame, int frame);


void inversebox(double a[3], double b[3], double c[3],
                double invbox[3][3]);
static void wrap_in_box(double a[3], double b[3], double c[3],
                        double invbox[3][3],
                        int natom, double x[MAXATOM], double y[MAXATOM],
                        double z[MAXATOM]);
void wrap_pairs(double *x, double *y, double *z,
                double a[3], double b[3], double c[3],
                int mol_list[MAXATOM][MAXATOM], int nmolecule,
                int bond_list[MAXATOM][MAXATOM], int natom, int j, int k,
                int wrapped[MAXATOM]);
int ele_index(char ea);
double atm_mass(char ea);
void read_bonds(double rbond[MAXELE][MAXELE], int *duration, int *xyz_copies, double *time_step,
                int *dump_structures);
double boxvol(double a[3], double b[3], double c[3]); // comments: simulate infinite environment with box periodicity
int retrieve_molecule(int atom_id, int mol_list[MAXATOM][MAXATOM]);


//int is_in_molecule(int atom_id, int mol_id, int mol_list[MAXATOM][MAXATOM]);

double xold[MAXHIST][MAXATOM], yold[MAXHIST][MAXATOM], zold[MAXHIST][MAXATOM];
int bond_duration = 8;  /* Determines how many frames a bond must exist. Default: 8 */
double rbond[MAXELE][MAXELE]; // comments: bond length data
int bond_list_old[MAXHIST][MAXATOM][MAXATOM];
int bond_curr[MAXATOM][MAXATOM]; /* whether atoms within bond distance for current timeframe only, no duration info */
int bond_prev[MAXATOM][MAXATOM];

// Define matrix of bond changes (breaking/forming) and count
int bond_changes[MAXATOM][MAXATOM][3];
int bond_changes_count;
// Define matrix of reactant/product molecules and counts
int involved[MAXATOM][2][MAXATOM];
int involved_count[MAXATOM][2];
int involved_reac_count;
// Define previous mol_list (array of atoms for each molecule)
int mol_list_prev[MAXATOM][MAXATOM];
int bonding_count = 0;
// Define matrix of atoms for each molecule in each reaction
int local_atoms[MAXATOM][2];
// Define array of molecule lifetimes and count
int lifetime[3000][MAXATOM];
int lifetime_count[3000];
// Define name size and array of molecule name
#define NAMESIZE 50
char lifetime_names[3000][MAXATOM][NAMESIZE];
// Define total number of molecules
int total_molecules = 0;
// Whether to count molecules
int count_molecules = 0;

int main(int argc, char** argv)
{
    char buf[BUFSIZE];
    double x[MAXATOM], y[MAXATOM], z[MAXATOM], q[MAXATOM];
    double r2;
    double rcomx[MAXATOM], rcomy[MAXATOM], rcomz[MAXATOM];
    double a[3], b[3], c[3];
    double invbox[3][3];
    int nmolecule;
    int natom, i, j, k, type[MAXATOM], mol_list[MAXATOM][MAXATOM],
        bond_list[MAXATOM][MAXATOM];
    char element[MAXELE];
    FILE *fout;
    FILE *fgen; // The gen (coordinate) file.
    FILE *fchr; // The chr (charge) file.
    int read_charge = 0; // Whether charges should be read.
    int natom1; // Number of atoms in the chr file.
    double dx, dy, dz;
    int xyz_copies = 0; // Print this many copies in each direction for the xyz file.
    int frame = 1;
    int ix, iy, iz;
    int printed_atom[MAXATOM]; // Whether an atom was printed out.
    double vbox; // The volume of the box.
    double time_step; // The time step between saved frames.
    int dump_structures = 0; // Whether to dump structures.
    int do_wrap_atoms = 1; // Whether to wrap atoms. This is required for correct molecule identification.
    int print_xyz = 0; // Whether to create molanal.xyz file.
    int bond_count[MAXELE][MAXELE];
    int atom_config[MAXATOM];
    int test = -1;
    int reaction_frame[MAXATOM][5];

    // Just defining all the variables which are going to be used in loops
	int ii = 0;
	int ll = 0;
	int kk = 0;
	int aa = 0;
	int bb = 0;
	int jj = 0;

	// Initializing atom_config
    for(ll= 0; ll<MAXATOM; ll++){
        atom_config[ll] = -1;
    }

    // Get data from bonds.dat: threshold of the bond distances between atoms (rbond), the bond duration to be considered as bonded (bond_duration),
    // the time step that a frame corresponds to (time_step), and data about dumping or not output files (xyz_copies and dump_structures)
    read_bonds(rbond, &bond_duration, &xyz_copies, &time_step, &dump_structures);

    // Output files that I never used
    if ( dump_structures ) {
        system("rm -rf molecules");
        system("mkdir molecules");
    }

    // Making rbond symmetric
    for ( i = 0; i < MAXELE; i++ ) {
        for ( j = i + 1; j < MAXELE; j++ ) {
            rbond[i][j] = rbond[j][i];
        }
    }


    if ( bond_duration >= MAXHIST ) {
        printf("Error: bond duration is too big.  Recompile with a bigger MAXHIST\n");
        exit(1);
    }

    // Initialize bond_prev to be NOT bonded
    // Initialize previous molecule list to -1
    for ( i = 0; i < MAXATOM; i++ ) {
        for ( j = 0; j < MAXATOM; j++ ) {
            bond_prev[i][j] = 0;
            mol_list_prev[i][j] = -1;
        }
    }

    //Initialize the lifetimes to -1 and the names
    for (i = 0; i < 3000; i++) {
        for (j = 0; j < MAXATOM; j++) {
            lifetime[i][j] = -1;
            memset(lifetime_names[i][j], 0, NAMESIZE);
        }
        lifetime_count[i] = 0;
    }

    // Initialize the old bond list
    for ( i = 0; i < MAXHIST; i++ ) {
        for ( j = 0; j < MAXATOM; j++ ) {
            for ( k = 0; k < MAXATOM; k++ ) {
                bond_list_old[i][j][k] = 0;
            }
        }
    }

    //Initialize the old coordinates of the atoms
    for ( i = 0; i < MAXATOM; i++ ) {
        q[i] = 0.0;
        for ( j = 0; j < bond_duration; j++ ) {
            xold[j][i] = NOTSET;
            yold[j][i] = NOTSET;
            zold[j][i] = NOTSET;
        }
    }

    if ( argc >= 2 ) {
        fgen = fopen(argv[1],"r");
            if ( fgen == NULL ) {
                printf("Error: could not open file %s\n", argv[1]);
                exit(1);
            }
    } else {
        printf("Error: The input file must be specified\n");
        exit(1);
    }

    if ( argc == 3 ) {
        read_charge = 1;
        fchr = fopen(argv[2],"r");
        if ( fchr == NULL ) {
            printf("Error: could not open file %s\n", argv[2]);
            exit(1);
        }
    }
    if (print_xyz) {
        fout = fopen("molanal.xyz", "w");
    }

    // Get the elements; THIS IS CURRENTLY A HACK
    element[0] = 'C';
    element[1] = 'H';

    // Reading through the .atom file
    while ( !feof(fgen) || ferror(fgen) ) {
        // Initializing the atoms coordinates and type
        for ( j = 0; j < MAXATOM; j++ ) {
            x[j] = 0.0;
            y[j] = 0.0;
            z[j] = 0.0;
            q[j] = 0.0;
            rcomx[j] = -100.0;
            rcomy[j] = -100.0;
            rcomz[j] = -100.0;
            type[j] = 0;
        }

       // Clear matrix of reactant/product molecules and counts and local_atoms
        for ( i = 0; i < MAXATOM; i++ ) {
            for ( j = 0; j < 2; j++ ) {
                for ( k = 0; k < MAXATOM; k++ ) {
                    involved[i][j][k] = -1;
                }
                involved_count[i][j] = 0;
                local_atoms[i][j] = -1;
            }
        }
        involved_reac_count = 0;


        // Reset bond changes matrix and count
        for ( i = 0; i < MAXATOM; i++ ) {
            for ( j = 0; j < MAXATOM; j++ ) {
                for ( k = 0; k < 3; k++ ) {
                    bond_changes[i][j][k] = -1;
                }
            }
   		}
        bond_changes_count = 0;

        if ( argc == 3 ) {
            for ( i = 0; i < 4; i++ ) {
                fgets(buf,BUFSIZ,fchr);
                if ( feof(fchr) ) {
                    printf("End of charge file reached\n");
                    exit(0);
                }
                if ( buf[0] != '#' ) {
                    printf("Bad format in charge file.\n");
                    exit(1);
                }
            }
            fscanf(fchr,"#   %d\n", &natom1);
            for ( i = 0; i < natom1; i++ ) {
                fgets(buf,BUFSIZ,fchr);
                if ( sscanf(buf,"%d  %lf", &j, &q[i]) != 2 ) {
                    printf("Error: could not read in a charge\n");
                    exit(1);
                }
            }
            if ( feof(fchr) || ferror(fchr) ) {
                printf("Error: end of CHR file found\n");
                exit(1);
            }
        }

        // Ignore Header
        fgets(buf,BUFSIZ,fgen); // ITEM: TIMESTEP
        fgets(buf,BUFSIZ,fgen); // timestep
        fgets(buf,BUFSIZ,fgen); // ITEM: NUMBER OF ATOMS

        // Find and store the number of atoms
        fgets(buf,BUFSIZ,fgen); // natoms
        if ( feof(fgen) || ferror(fgen) ) {
            break;
        }
        sscanf(buf,"%d", &natom);
        if ( read_charge && natom != natom1 ) {
            printf("Error: number of atoms in the chr and gen files are different\n");
            exit(1);
        }
        if ( natom > MAXATOM ) {
            printf("Error: the number of atoms is too large.  Increase MAXATOM and recompile\n");
        }

        // Parse the box Bounds
        fgets(buf,BUFSIZ,fgen); // ITEM: BOX BOUNDS

        fgets(buf,BUFSIZ,fgen);
        sscanf(buf,"%lf %lf", &a[1], &a[0]);
        a[0] = a[0] - a[1];
        a[1] = 0.0;
        a[2] = 0.0;
        fgets(buf,BUFSIZ,fgen);
        sscanf(buf,"%lf %lf", &b[0], &b[1]);
        b[1] = b[1] - b[0];
        b[0] = 0.0;
        b[2] = 0.0;
        fgets(buf,BUFSIZ,fgen);
        sscanf(buf,"%lf %lf", &c[0], &c[2]);
        c[2] = c[2] - c[0];
        c[0] = 0.0;
        c[1] = 0.0;
        if ( a[1] > 0.0 || b[0] > 0.0 || c[2] <= 0.0 ) {
            printf("Bad cell size read in\n");
            exit(1);
        }

        fgets(buf,BUFSIZ,fgen); // ITEM: ATOMS id type xs ys zs

        /** Parse the input line **/
        for ( j = 0; j < natom; j++ ) {
            fgets(buf,BUFSIZ,fgen);
            sscanf(buf,"%d",&i);                                                                                                                                                                         // associate line with particular atom
            sscanf(buf,"%d %d %lf %lf %lf", &i, &type[i-1], &x[i-1], &y[i-1], &z[i-1]);

            /** Un-scale the atom coordinates **/
            x[i-1] = x[i-1]*a[0];
            y[i-1] = y[i-1]*b[1];
            z[i-1] = z[i-1]*c[2];
        }

        /** Wrap all coordinates back into the primitive unit cell. */
        inversebox(a,b,c,invbox);
        wrap_in_box(a,b,c,invbox,natom,x,y,z);
        vbox = boxvol(a,b,c);

        // Initialize mol_list (1st index represents a molecule and the 2nd index are the indexes of the atoms in this molecule),
        // bond_list (bond_list[i][j] = bond_list[j][i] = 1 if i and j are bonded, 0 otherwise), bond_curr is kind of the same (??)
        for ( j = 0; j < MAXATOM; j++ ) {
            for ( i = 0; i < MAXATOM; i++ ) {
                mol_list[i][j] = -1;
                bond_list[i][j] = 0;
                bond_curr[i][j] = 0;
           }
        }

        nmolecule = 0;

        /* Compute bond list by comparing the size of the bonds to the associated bond threshold */
        for (j=0; j<natom-1; j++) {
            for (i=j+1; i<natom; i++) {
                r2 = wrap_atom(j,i,x,y,z,a,b,c,MAXWRAP,0);
                if (r2<r2bond(element[type[j]-1], element[type[i]-1])) {
                    bond_curr[j][i] = 1;
                    bond_curr[i][j] = 1; // comments: bond_curr is based on only bond length criteria
                }
            }
        }

        /** Each atom is assigned to a molecule **/
        find_molecules(x,y,z,a,b,c,mol_list,natom,&nmolecule,
               type,element,bond_list);

        // Don't output the frames before bond_duration because we don't consider there is any bonds yet
        if ( frame > bond_duration ) {
            printf("Beginning frame %d\n", frame - bond_duration );
        }
        total_molecules += nmolecule;
        /** Wrap each atom to make whole molecules **/
        /** wrap_atoms(x, y, z, a, b, c, bond_list, natom) ; **/
        wrap_molecules(x, y, z, a, b, c, mol_list, nmolecule,bond_list, natom);

        /** Wrap the center of mass of each molecule to be inside the box**/
        wrap_com(x,y,z,a,b,c,mol_list,nmolecule, natom,
     			         rcomx, rcomy, rcomz,type,element);


        for(ii=0; ii<MAXATOM; ii++){
            for(jj = 0; jj<5;jj++){
                reaction_frame[ii][jj] = -1;
            }
        }

        //
        int counter = 0;
        for (ii=0;ii<natom;ii++){
            // Consider one atom by one
            int atom = ii;
            int molecule = retrieve_molecule(ii, mol_list);
            //Initialize bond_count which is going to ouput the configuration (or features) of the atom
            for (ll = 0; ll< MAXELE; ll++ ) {
                for (kk = 0; kk < MAXELE; kk++ ) {
                    bond_count[ll][kk] = 0;
                }
            }

            // Count the number of bonds of each type for the atom in question
            for (ll = 0; mol_list[molecule][ll] >= 0; ll++ ) {
                int jj = mol_list[molecule][ll];
                if ( bond_list[atom][jj] == 1) {
                    if ( type[atom] > type[jj] ) {
                        bond_count[type[atom]-1][type[jj]-1]++;
                    } else {
                        bond_count[type[jj]-1][type[atom]-1]++;
                    }
                }
            }
            // This is creating the configuration of the atom ( 1004, 2001....)
            int local_counts[3] = {0, 0, 0}; // (C-C), (H-H), (H-C)

            //Store the number of bonds of each type in local_counts for this atom
            for (ll = 0; ll < MAXELE; ll++ ) {
                for (kk = 0; kk < MAXELE; kk++ ) {
                    if ( bond_count[kk][ll] > 0 ) {
                        if (kk == 0 && ll == 0) local_counts[0] += bond_count[kk][ll];
                        if (kk == 1 && ll == 1) local_counts[1] += bond_count[kk][ll];
                        if (kk == 1 && ll == 0) local_counts[2] += bond_count[kk][ll];
                    }
                }
            }

            // test is the configurations of the atom.
            test = type[ii]*1000 + local_counts[0]*100 + local_counts[1]*10 + local_counts[2];
            int tester = 0;

            for (ll =0; ll<MAXATOM; ll++){
                // Test if a reaction has happened to the atom by seeing if its atom configuration changes or if it is involved in a bond change.
                if (test != atom_config[ii] || bond_changes[molecule][ll][1] == ii ||bond_changes[molecule][ll][2] == ii ){
                    // You just to see one change in the atom to store it.
                    if (tester == 0){
                        // Store the info of this atom in reaction_frame: the index of the atom, its old and new atom configuration,
                        // its old and new molecule.
                        reaction_frame[counter][0] = ii;
                        reaction_frame[counter][1] = atom_config[ii];
                        reaction_frame[counter][2] = test;
                        reaction_frame[counter][3] = retrieve_molecule(ii, mol_list_prev);
                        reaction_frame[counter][4] = molecule;
                        // Print the atoms bonded to this atom to be able to initialize the system
                        if (frame==1){
                            printf("ID: %d| %d \n", atom, test);
                            printf("Bond: %d is bonded to ", atom);
                            for (jj = 0; jj<MAXATOM; jj++){
                                if (bond_list[jj][ii] == 1){
                                    printf("%d ", jj);
                                }
                            }
                            printf("\n");
                        }
                        //Store the atom configuration in atom_config
                        atom_config[ii] = test;
                        // Counter counts the number of atoms which are modified in this frame
                        counter++;
                        tester++;
                    }
                }
            }
        }

        // This is where you're going to print the reactions. You don't do it at the first frame.
        if(frame>1){
            // Store the old and new atoms configuration and their associated index in sort_help.
            int sort_featfeat[counter][5];
            int sort_help[counter][3];
            int indexes[counter];
            for(jj=0; jj<counter; jj++){
                sort_help[jj][0] = reaction_frame[jj][1];
                sort_help[jj][1] = jj;
                sort_help[jj][2] = reaction_frame[jj][2];
            }

            // Sort the atoms which have reacted by order of configuration (1001<1004<2001). It is important
            // to keep this order to avoid counting several times the reactions for example.
            for(jj=0; jj<counter; jj++){
                for(kk=0; kk<counter-1; kk++){
                    if(sort_help[kk][0]>sort_help[kk+1][0]){
                        int temp = sort_help[kk+1][0];
                        sort_help[kk+1][0] = sort_help[kk][0];
                        sort_help[kk][0] = temp;
                        temp = sort_help[kk+1][1];
                        sort_help[kk+1][1] = sort_help[kk][1];
                        sort_help[kk][1] = temp;
                        temp = sort_help[kk+1][2];
                        sort_help[kk+1][2] = sort_help[kk][2];
                        sort_help[kk][2] = temp;
                    }
                    else if(sort_help[kk][0]==sort_help[kk+1][0]) {
                        if (sort_help[kk][2]>sort_help[kk+1][2]){
                            int temp = sort_help[kk+1][0];
                            sort_help[kk+1][0] = sort_help[kk][0];
                            sort_help[kk][0] = temp;
                            temp = sort_help[kk+1][1];
                            sort_help[kk+1][1] = sort_help[kk][1];
                            sort_help[kk][1] = temp;
                            temp = sort_help[kk+1][2];
                            sort_help[kk+1][2] = sort_help[kk][2];
                            sort_help[kk][2] = temp;
                        }
                    }
                }
            }

            // Store all the info from reaction_frame in the right order in sort_featfeat
            for(jj=0;jj<counter;jj++){
                for(kk=0;kk<5;kk++){
                    sort_featfeat[jj][kk] = reaction_frame[sort_help[jj][1]][kk];
                }
            }

            for(jj=0;jj<counter;jj++){
                for(kk=0;kk<5;kk++){
                    reaction_frame[jj][kk] = sort_featfeat[jj][kk];
                }
            }

            // Initialize molecule_involved
            int molecule_involved[counter][2];
            for(kk=0; kk<counter;kk++){
                for(ll=0; ll<2; ll++){
                    molecule_involved[kk][ll] = -1;
                }
            }

            //Initialize bond_in_reactant and bond_in_product
            int bond_in_reactant[counter][counter];
            int bond_in_product[counter][counter];
            for (kk = 0; kk<counter; kk++){
                for (ll = 0; ll<counter; ll++){
                    bond_in_reactant[kk][ll] = -1;
                    bond_in_product[kk][ll] = -1;
                }
            }

            // A bond is either broken or created. If the bond is broken, the bond was present in the reactants, otherwise
            // the bond is present in the products.
            for (kk = 0; kk<counter; kk++){
                for (ll = kk +1; ll<counter; ll++){
                    if (bond_prev[reaction_frame[kk][0]][reaction_frame[ll][0]] == 1){
                        bond_in_reactant[kk][ll] = 1;
                        bond_in_reactant[ll][kk] = 1;
                    }
                    if (bond_list[reaction_frame[kk][0]][reaction_frame[ll][0]] == 1){
                        bond_in_product[kk][ll] = 1;
                        bond_in_product[ll][kk] = 1;
                    }
                }
            }

            // Print the list of all the atoms involved in reactions in this frame.
            if (counter>0){
                printf("Atoms involved: ");
                for(kk = 0; kk<counter; kk++){
                    printf("%d ", reaction_frame[kk][0]);
                }
                printf("\n");
            }

            // Store index and config of each involved atoms in local_reaction_frame.
            int local_reaction[2][5];
            int local_reaction_frame[counter][2];
            int type_of_bond;
            for (ll = 0; ll<counter; ll++){
                local_reaction_frame[ll][0] = reaction_frame[ll][0];
                local_reaction_frame[ll][1] = reaction_frame[ll][1];
            }

            // Print the bonds which change and the reaction written with the local features.
            for(ll = 0; ll<MAXATOM; ll++){
                for (jj = ll+1; jj<MAXATOM; jj++){
                    // A reaction has happened if the bond in the previous frame is not the same as the bond in this frame.
                    if (bond_list[ll][jj] != bond_prev[ll][jj]){
                        int local_counter = 0;
                        for (kk = 0; kk<counter; kk++){
                            // Store in local_reaction the info about the reactants.
                            if (ll == local_reaction_frame[kk][0]){
                                local_reaction[local_counter][0] = local_reaction_frame[kk][0];
                                local_reaction[local_counter][1] = local_reaction_frame[kk][1];
                                local_reaction[local_counter][2] = type[ll];
                                local_reaction[local_counter][4] = kk;
                                local_counter++;
                            }
                            if (jj == local_reaction_frame[kk][0]){
                                local_reaction[local_counter][0] = local_reaction_frame[kk][0];
                                local_reaction[local_counter][1] = local_reaction_frame[kk][1];
                                local_reaction[local_counter][2] = type[jj];
                                local_reaction[local_counter][4] = kk;
                                local_counter++;
                            }
                        }
                        // Calculate the new configuration of the atoms due to this reaction (which doesn't have to be the
                        // configuration at the end of the frame)
                        if (type[jj]+type[ll] == 4){
                            type_of_bond = 10;
                        }
                        else if(type[jj]+type[ll] == 3){
                            type_of_bond = 1;
                        } else{
                            type_of_bond = 100;
                        }
                        if(bond_list[ll][jj] == 1){
                            local_reaction[0][3] = local_reaction[0][1] + type_of_bond;
                            local_reaction[1][3] = local_reaction[1][1] + type_of_bond;
                        }
                        if(bond_list[ll][jj] == 0){
                            local_reaction[0][3] = local_reaction[0][1] - type_of_bond;
                            local_reaction[1][3] = local_reaction[1][1] - type_of_bond;
                        }
                        // Print the change in bond and the reaction with local features. The latter has a specific order:
                        // the first reactant of the reaction is the one with the lowest configuration and the order is kept for the products.
                        if (local_reaction[0][1]<local_reaction[1][1]){
                            printf("Change bond: %d and %d from %d to %d\n", local_reaction[0][0],local_reaction[1][0],bond_prev[ll][jj], bond_list[ll][jj]);
                            printf("Local reaction: %d + %d => %d + %d \n", local_reaction[0][1], local_reaction[1][1], local_reaction[0][3], local_reaction[1][3]);
                        }
                        else if (local_reaction[0][1]>local_reaction[1][1]){
                            printf("Change bond: %d and %d from %d to %d\n", local_reaction[1][0], local_reaction[0][0],bond_prev[ll][jj], bond_list[ll][jj]);
                            printf("Local reaction: %d + %d => %d + %d \n", local_reaction[1][1], local_reaction[0][1], local_reaction[1][3], local_reaction[0][3]);
                        }
                        else{
                            if (local_reaction[0][3]<local_reaction[1][3]){
                                printf("Change bond: %d and %d from %d to %d\n", local_reaction[0][0], local_reaction[1][0],bond_prev[ll][jj], bond_list[ll][jj]);
                                printf("Local reaction: %d + %d => %d + %d \n", local_reaction[0][1], local_reaction[1][1], local_reaction[0][3], local_reaction[1][3]);
                            }
                            else{
                                printf("Change bond: %d and %d from %d to %d\n", local_reaction[1][0], local_reaction[0][0],bond_prev[ll][jj], bond_list[ll][jj]);
                                printf("Local reaction: %d + %d => %d + %d \n", local_reaction[1][1], local_reaction[0][1], local_reaction[1][3], local_reaction[0][3]);
                            }
                        }
                        local_reaction_frame[local_reaction[0][4]][1] = local_reaction[0][3];
                        local_reaction_frame[local_reaction[1][4]][1] = local_reaction[1][3];
                        if(bond_list[ll][jj] == 1){
                            printf("Local bond reaction: %d + %d => %d + %d \n", 0, 0, 2, 1);
                        }
                        if(bond_list[ll][jj] == 0){
                            printf("Local bond reaction: %d + %d => %d + %d \n", 2, 1, 0, 0);
                        }
                    }

                }
            }

            int m_list[MAXATOM][MAXATOM];
            int b_list[MAXATOM][MAXATOM];

            if (counter !=0){
                // Print the reaction with all the atoms for this frame. I don't use it anymore, it can probably be deleted.
                // (only the 'Reaction with atoms' part and the 'Reaction with bond' part can be deleted). The rest is useful.
                printf("Reaction with atoms: ");
                    for (kk = 0; kk<counter;kk++){
                        printf("%d ", reaction_frame[kk][1]);
                        if(kk!=counter-1){
                            printf("+ ");
                        }
                        else{
                            printf("=> ");
                        }
                    }

                    for (kk = 0; kk<counter;kk++){
                        printf("%d ", reaction_frame[kk][2]);
                        if(kk!=counter-1){
                            printf("+ ");
                        }
                        else{
                            printf("\n");
                        }
                    }
                printf("Reaction with bond: ");
                for (kk = 0; kk<counter; kk++){
                    int count_bond = 0;
                    for(ll = 0; ll<counter; ll++){
                        if (bond_in_reactant[kk][ll]>0){
                            printf("%d ", ll + 1);
                            count_bond++;
                        }
                    }
                    if (count_bond ==0){
                        printf("0 ");
                    }
                    if (kk != counter-1){
                        printf("+ ");
                    } else{
                        printf("=> ");
                    }
                }

                for (kk = 0; kk<counter; kk++){
                    int count_bond = 0;
                    for(ll = 0; ll<counter; ll++){
                        if (bond_in_product[kk][ll]>0){
                            printf("%d ", ll + 1);
                            count_bond++;
                        }
                    }
                    if (count_bond ==0){
                        printf("0 ");
                    }
                    if (kk != counter-1){
                        printf("+ ");
                    } else{
                        printf("\n");
                    }
                }


                //This is to print the reaction involving the molecules for the whole frame. It is kept because it is
                // what is used for the 'MD' graph in the resutls.
                printf("Reaction with molecules: ");
                int mol_react = 0;
                int mol_prod = 0;
                // Store what are the molecules which are modified during this frame.
                for(kk = 0; kk<counter; kk++){
                    if(mol_react==0){
                        molecule_involved[mol_react][0] = reaction_frame[kk][3];
                        mol_react++;
                    } else {
                        int tester_1 = 0;
                        for (ll=0; ll<mol_react;ll++){
                            if(molecule_involved[ll][0] == reaction_frame[kk][3]){
                                tester_1++;
                            }
                        }
                        if(tester_1 == 0){
                            molecule_involved[mol_react][0] =reaction_frame[kk][3];
                            mol_react++;
                        }
                    }

                    if(mol_prod==0){
                        molecule_involved[mol_prod][1] = reaction_frame[kk][4];
                        mol_prod++;
                    } else {
                        int tester_2 = 0;
                        for (ll=0; ll<mol_prod;ll++){
                            if(molecule_involved[ll][1] == reaction_frame[kk][4]){
                                tester_2++;
                            }
                        }
                        if(tester_2 == 0){
                            molecule_involved[mol_prod][1] = reaction_frame[kk][4];
                            mol_prod++;
                        }
                    }
                }

                // Print the reaction, by first writing the molecule in the normal convention (C1 H4 4(H-C), H2 1(H-H)...)
                int num_mol;
                int ele_counts[MAXATOM][MAXELE];
                int bond_counts[MAXATOM][MAXELE][MAXELE];
                for(ii = 0; ii<2; ii++){
                    if(ii == 0){
                        num_mol = mol_react;
                    } else {
                        num_mol = mol_prod;
                    }

                    if (ii == 1) {
                        memcpy(m_list, mol_list, sizeof m_list);
                        memcpy(b_list, bond_list, sizeof b_list);
                    } else {
                        memcpy(m_list, mol_list_prev, sizeof m_list);
                        memcpy(b_list, bond_prev, sizeof b_list);
                    }

                    for(ll=0; ll<num_mol; ll++){
                        int molecule = molecule_involved[ll][ii];

                        for (aa = 0; aa < MAXELE; aa++ ) {
                            ele_counts[ll][aa] = 0;
                        }

                        for (aa = 0; m_list[molecule][aa] >= 0; aa++ ) {
                            ele_counts[ll][type[m_list[molecule][aa]]-1]++;
                        }

                        // Reset bond counts
                        for ( aa = 0; aa < MAXELE; aa++ ) {
                            for (bb = 0; bb < MAXELE; bb++ ) {
                                bond_counts[ll][aa][bb] = 0;
                            }
                        }

                        // Calculate bond counts
                        for ( aa = 0; m_list[molecule][aa] >= 0; aa++ ) {
                            int kk = m_list[molecule][aa];
                            for (bb = 0; m_list[molecule][bb] >= 0; bb++ ) {
                                int jj = m_list[molecule][bb];
                                if ( b_list[kk][jj] == 1 && jj > kk ) {
                                    if ( type[kk] > type[jj] ) {
                                        bond_counts[ll][type[kk]-1][type[jj]-1]++;
                                    } else {
                                        bond_counts[ll][type[jj]-1][type[kk]-1]++;
                                    }
                                }
                            }
                        }

                        for(kk = 0; kk<MAXELE; kk++){
                            if (ele_counts[ll][kk]!=0){
                                printf("%c%d ", element[kk], ele_counts[ll][kk]);
                            }
                        }
                        for(kk = 0; kk<MAXELE;kk++){
                            for(jj = 0; jj<MAXELE; jj++){
                                if (bond_counts[ll][kk][jj]!=0){
                                    printf("%d(%c-%c) ", bond_counts[ll][kk][jj], element[kk], element[jj]);
                                }
                            }
                        }

                        if(ll!=num_mol - 1){
                            printf("+ ");
                        } else {
                            if(ii == 0){
                                printf("=> ");
                            } else {
                                printf("\n");
                            }
                        }
                    }
                }
            }
        }

        // For the first frame, print out what are the molecules so you can follow them more easily after.
        int ele_counts[MAXATOM][MAXELE];
        int bond_counts[MAXATOM][MAXELE][MAXELE];
        if (frame==1){
            for(ll=0; ll<nmolecule; ll++){
                int molecule = ll;

                for (aa = 0; aa < MAXELE; aa++ ) {
                    ele_counts[ll][aa] = 0;
        		}

                for (aa = 0; mol_list[molecule][aa] >= 0; aa++ ) {
                    ele_counts[ll][type[mol_list[molecule][aa]]-1]++;
                }

                // Reset bond counts
                for (aa = 0; aa < MAXELE; aa++ ) {
                    for (bb = 0; bb < MAXELE; bb++ ) {
                        bond_counts[ll][aa][bb] = 0;
                    }
        		}

                // Calculate bond counts
        		for (aa = 0; mol_list[molecule][aa] >= 0; aa++ ) {
                    int kk = mol_list[molecule][aa];
                    for (bb = 0; mol_list[molecule][bb] >= 0; bb++ ) {
                        int jj = mol_list[molecule][bb];
                        if ( bond_list[kk][jj] == 1 && jj > kk ) {
                            if ( type[kk] > type[jj] ) {
                                bond_counts[ll][type[kk]-1][type[jj]-1]++;
                            } else {
                                bond_counts[ll][type[jj]-1][type[kk]-1]++;
                            }
                        }
                    }
        		}
                printf("Molecule: ");
        		for(kk = 0; kk<MAXELE; kk++){
        		    if (ele_counts[ll][kk]!=0){
        		        printf("%c%d ", element[kk], ele_counts[ll][kk]);
        		    }
        		}
        		for(kk = 0; kk<MAXELE;kk++){
        		    for(jj = 0; jj<MAXELE; jj++){
        		        if (bond_counts[ll][kk][jj]!=0){
        		            printf("%d(%c-%c) ", bond_counts[ll][kk][jj], element[kk], element[jj]);
        		        }
        		    }
        		}
        		printf("\n");
            }
        }


        //Not used because I don't output an xyz file usually.
        if ( frame > bond_duration ) {
            for ( j = 0; j < natom; j++ ) {
                printed_atom[j] = 0;
            }
            bonding_count = 0;

            if (print_xyz) {
                fprintf(fout,"%d\n\n", natom * (2*xyz_copies+1) * (2*xyz_copies+1) * (2*xyz_copies+1));
                for ( ix = -xyz_copies; ix <= xyz_copies; ix++ ) {
                    for ( iy = -xyz_copies; iy <= xyz_copies; iy++ ) {
                        for ( iz = -xyz_copies; iz <= xyz_copies; iz++ ) {
                            dx = ix * a[0] + iy * b[0] + iz * c[0];
                            dy = ix * a[1] + iy * b[1] + iz * c[1];
                            dz = ix * a[2] + iy * b[2] + iz * c[2];
                            for ( j = 0; j < natom; j++ ) {
                                fprintf(fout,"%c %f %f %f\n", element[type[j]-1],
                                dx + x[j], dy + y[j], dz + z[j]);
                            }
                        }
                    }
                }
            }
            printf("End frame %d\n", frame- bond_duration);
        }

        if ( feof(fgen) || ferror(fgen) ) {
                break;
        }

        // Update the positions of atoms.
        for ( j = 0; j < natom; j++ ) {
            for ( i = 0; i < bond_duration - 1; i++ ) {
                xold[i][j] = xold[i+1][j];
                yold[i][j] = yold[i+1][j];
                zold[i][j] = zold[i+1][j];
            }
            if ( bond_duration > 0 ) {
                xold[bond_duration-1][j] = x[j];
                yold[bond_duration-1][j] = y[j];
                zold[bond_duration-1][j] = z[j];
            }
        }
        /* Update bond_list history */
        for (j=0; j<natom-1; j++) {
            for (k=j+1; k<natom; k++) {
                for (i=0; i<bond_duration-1; i++) {
                    bond_list_old[i][j][k] = bond_list_old[i+1][j][k];
                    bond_list_old[i][k][j] = bond_list_old[i+1][k][j];
                }
                if (bond_duration>0) {
                    bond_list_old[bond_duration-1][j][k] = bond_curr[j][k]; // bond_list_old is defined only by bond length (from bond_curr)
                    bond_list_old[bond_duration-1][k][j] = bond_curr[k][j];
                }
            }
        }

        // bond_prev is the true bonded flag (bond length & duration)
        memcpy(bond_prev, bond_list, sizeof bond_prev);
        // Set mol_list_prev to the current mol_list before looping
        memcpy(mol_list_prev, mol_list, sizeof mol_list_prev);

        frame++;
    }
    if (print_xyz) {
        fclose(fout);
    }

    if (count_molecules) {
        printf("\n ----------------------------------------\n");
        printf("Total Number of Molecules: %d\n", total_molecules);
        for (i = 0; i < 3000; i++) {
            printf("Frame %d\n", i);
            for (j = 0; lifetime[i][j] >= 0; j++) {
                printf("%s: %d\n", lifetime_names[i][j], lifetime[i][j]);
            }
        }
    }

    return(0);
}

void find_molecules(double *x, double *y, double *z,
                    double a[3], double b[3], double c[3],
                    int mol_list[MAXATOM][MAXATOM],
                    int natom, int *nmolecule,
                    int *type, char *element, int bond_list[MAXATOM][MAXATOM])
{
    int j, k, l;
    int sort_list[MAXATOM];
    int used[MAXATOM];
    int best_idx;
    int best;

    for ( j = 0; j < natom; j++ ) {
        if ( not_in_molecule(j,mol_list,*nmolecule) ) {
            mol_list[*nmolecule][0] = j;
            find_molecule(j,x,y,z,a,b,c,mol_list,natom,*nmolecule,
                  type, element,bond_list);
            #if (0)
            {
                int k;
                for ( k = 0; mol_list[*nmolecule][k] >= 0; k++ ) {
                    printf("mol:%d %d %d\n", *nmolecule, k, mol_list[*nmolecule][k]);
                }
            }
            #endif
            (*nmolecule)++;
            bond_changes_count = 0;
        }
    }
}



void wrap_atoms(double *x, double *y, double *z,
                double a[3], double b[3], double c[3],
                int bond_list[MAXATOM][MAXATOM],
                int natom)
/* Wrap bonded atoms so that they are as close as possible to each other. */
{
    int j, k;
    for ( j = 0; j < natom; j++ ) {
        for ( k = j+1; k < natom; k++ ) {
            if ( bond_list[j][k] == 1 ) {
                wrap_atom(j, k, x,y,z,a,b,c,MAXWRAP,1);
            }
        }
    }
}


int not_in_molecule(int isearch, int mol_list[MAXATOM][MAXATOM],
                    int nmolecule)
/** Return 1 if atom isearch is not in the given molecule, 0
   otherwise **/
{
    int j, k;
    for ( j = 0; j < nmolecule; j++ ) {
        for ( k = 0; mol_list[j][k] >= 0; k++ ) {
            if ( mol_list[j][k] == isearch ) {
                return ( 0 );
            }
        }
    }
    return (1);
}

void find_molecule(int j,
                   double *x, double *y, double *z,
                   double a[3],  double b[3], double c[3],
                   int mol_list[MAXATOM][MAXATOM],
                   int natom,
                   int nmolecule,
                   int *type,
                   char *element,
                   int bond_list[MAXATOM][MAXATOM])
/** Recursively search for all atoms bonded to j, then all atoms
    bonded to atoms bonded to j, etc. **/
{
    int k, l;
    for ( k = 0; k < natom; k++ ) {
        if ( is_bonded(j,k,x,y,z,a,b,c,type,element,mol_list,nmolecule) ) {
            bond_list[j][k] = 1; // comments: true bonded flag set here
            for ( l = 0; mol_list[nmolecule][l] >= 0; l++ ) {
                if ( mol_list[nmolecule][l] == k ) break;
            }
            if ( mol_list[nmolecule][l] < 0 ) {
                /** We found a new entry **/
                mol_list[nmolecule][l] = k;
                find_molecule(k,x,y,z,a,b,c,mol_list,natom,nmolecule,
                          type,element,bond_list);
            }
        }
    }
}


int is_bonded(int j, int k, double *x, double *y, double *z,
              double a[3], double b[3], double c[3],
              int *type, char *element, int mol_list[MAXATOM][MAXATOM], int nmolecule)
/** Return 1 if j and k are bonded, 0 if not **/
/** x, y, and z contain the atomic coordinates. **/
/** a, b, and c are the rectangular box dimensions **/
{
    // double r2 ;
    int i, h, t;
    if ( j == k ) return 0;
    #if (0)
		printf("\nTesting Atoms %d (%c) and %d (%c) for bonding: r = %f\n", j, element[type[j]-1], k, element[type[k]-1], sqrt(r2) );
		fflush(stdout);
		printf("i, j, k = %d %d %d\n", i, j, k);
		printf("bond_curr[j][k] = %d\n", bond_curr[j][k]);
		printf("bond_prev[j][k] = %d\n", bond_prev[j][k]);
		printf("xold[i][j] = %f\n", xold[i][j]);
		printf("bond_list_old[i][j][k] = %d\n", bond_list_old[i][j][k]);
    #endif

    /* Variable definitions
    - bond_curr = current connectivity graph (set by bond length data)
    - bond_prev = previous bonded graph (set from bond_list at end of main)
    - bond_list = current bonded graph (set by this function, is_bonded)
    - bond_list_old = past connectivity graph (set from bond_curr at end of main)
    - xold = previous x position? what is NOTSET? (set from x)
    */

    // The present passes the bonding criterion, now look at the past
    if ( bond_curr[j][k] ) { // Based on bond length criteria
        // If bond had not previously stabilized, then new bond requires bond duration to form
        if ( !bond_prev[j][k] ) { // Based on bond length criteria and connectivity of last 8 steps
            for ( i = 0; i < bond_duration; i++ ) {
                if ( xold[i][j] != NOTSET ) {
                    // Failing a bond test at any point in the past means that we are not bonded
                    if ( bond_list_old[i][j][k] == 0 ) { // comments: based on bond length criteria
                        return(0);
                    }
                }
            }
            // Bond formed
            if (j > k) {
                bond_changes[nmolecule][bond_changes_count][0] = 1;
                bond_changes[nmolecule][bond_changes_count][1] = j;
                bond_changes[nmolecule][bond_changes_count][2] = k;
                bond_changes_count++;
            }
        }
        return (1);
    }
    else {
        // If bond was previously stabilized, bond breaking requires duration to break
        if (bond_prev[j][k] ) {
            for ( i = 0; i < bond_duration; i++ ) {
                if ( xold[i][j] != NOTSET ) {
                    // Passing a bond test at any point in the past means that we are still bonded
                    if ( bond_list_old[i][j][k] ) {
                        return(1);
                    }
                }
            }
            // Bond broken
            if (j > k) {
                bond_changes[nmolecule][bond_changes_count][0] = 0;
                bond_changes[nmolecule][bond_changes_count][1] = j;
                bond_changes[nmolecule][bond_changes_count][2] = k;
                bond_changes_count++;
            }
        }
        return (0);
    }
}



double wrap_atom(int j, int k, double *x, double *y, double *z,
                 double a[3], double b[3], double c[3], int maxwrap,
                 int do_wrap)
/** Wrap atom k so to be as close as possible to atom j if do_wrap is
   non-zero.
   Returns the closest distance**2 found.  **/
{
    int l, m, n;
    double x1[3];
    double bestdist = 0.0, dist;
    int bestl = 0, bestm = 0, bestn = 0;

    bestdist = 1.0e20;

    for ( l = -maxwrap; l <= maxwrap; l++ ) {
        for ( m = -maxwrap; m <= maxwrap; m++ ) {
            for ( n = -maxwrap; n <= maxwrap; n++ ) {
                x1[0] = x[j] + l * a[0] + m * b[0] + n * c[0];
                x1[1] = y[j] + l * a[1] + m * b[1] + n * c[1];
                x1[2] = z[j] + l * a[2] + m * b[2] + n * c[2];

                dist = (x1[0]-x[k]) * (x1[0] - x[k]) +
                       (x1[1] - y[k]) * (x1[1] - y[k]) +
                       (x1[2] - z[k]) * (x1[2] - z[k]);
                if ( dist < bestdist ) {
                    bestdist = dist;
                    bestl = l;
                    bestm = m;
                    bestn = n;
                }
            }
        }
    }
    if ( do_wrap ) {
    #if (0)
        printf("bestdist = %f\n", bestdist);
        printf("bestl = %d, bestm = %d, bestn = %d\n", bestl, bestm, bestn);
        printf("dx = %f\n", bestl * a[0] + bestm * b[0] + bestn * c[0] );
        printf("dy = %f\n", bestl * a[1] + bestm * b[1] + bestn * c[1] );
        printf("dz = %f\n", bestl * a[2] + bestm * b[2] + bestn * c[2] );
    #endif
    x[k] -= bestl * a[0] + bestm * b[0] + bestn * c[0];
    y[k] -= bestl * a[1] + bestm * b[1] + bestn * c[1];
    z[k] -= bestl * a[2] + bestm * b[2] + bestn * c[2];
    }
    return bestdist;
}



void wrap_com(double *x, double *y, double *z,
              double a[3], double b[3], double c[3],
              int mol_list[MAXATOM][MAXATOM], int nmolecule, int natom,
              double rcomx[MAXATOM], double rcomy[MAXATOM],
              double rcomz[MAXATOM], int type[MAXATOM],
              char element[MAXELE])
/** Wrap each molecule so that it's center of mass (ignoring
   atomic weight) is inside the box **/
{
    double rx, ry, rz, dx, dy, dz, rx0, ry0, rz0, mtot;
    double mass;
    double x1[2], y1[2], z1[2];
    double mass_cell;

    int j, k, kk;

    rx0 = 0.0;
    ry0 = 0.0;
    rz0 = 0.0;
    mass_cell = 0.0;
    for ( j = 0; j < nmolecule; j++ ) {
        rx = 0.0;
        ry = 0.0;
        rz = 0.0;
        mtot = 0.0;
        for ( k = 0; mol_list[j][k] >= 0; k++ ) {
            kk = mol_list[j][k];
            mass = atm_mass(element[type[mol_list[j][k]]-1]);
            mtot += mass;
            rx += x[mol_list[j][k]] * mass;
            ry += y[mol_list[j][k]] * mass;
            rz += z[mol_list[j][k]] * mass;
        }
        /* rx, ry, rz is the molecular center of mass. */
        rx /= mtot;
        ry /= mtot;
        rz /= mtot;

        rcomx[j] = rx;
        rcomy[j] = ry;
        rcomz[j] = rz;

        x1[0] = rx;
        y1[0] = ry;
        z1[0] = rz;

        /** Calculate the center of the unit cell, with crystal coordinates
           (0.5, 0.5, 0.5) */
        x1[1] = 0.5 * a[0] + 0.5 * b[0] + 0.5 * c[0];
        y1[1] = 0.5 * a[1] + 0.5 * b[1] + 0.5 * c[1];
        z1[1] = 0.5 * a[2] + 0.5 * b[2] + 0.5 * c[2];

        /* Find the wrapped coordinates closest to the center of the unit cell. */
        wrap_atom(1,0,x1,y1,z1,a,b,c,MAXWRAP,1);

        dx = rx - x1[0];
        dy = ry - y1[0];
        dz = rz - z1[0];
        #if (1)
            rcomx[j] -= dx;
            rcomy[j] -= dy;
            rcomz[j] -= dz;

            for ( k = 0; mol_list[j][k] >= 0; k++ ) {
                x[mol_list[j][k]] -= dx;
                y[mol_list[j][k]] -= dy;
                z[mol_list[j][k]] -= dz;
            }
            for ( k = 0; mol_list[j][k] >= 0; k++ ) {
                mass = atm_mass(element[type[mol_list[j][k]]-1]);
                rx0 += x[mol_list[j][k]] * mass;
                ry0 += y[mol_list[j][k]] * mass;
                rz0 += z[mol_list[j][k]] * mass;
                mass_cell += mass;
            }
        #endif
    }
    dx = rx0 / mass_cell;
    dy = ry0 / mass_cell;
    dz = rz0 / mass_cell;
    for ( j = 0; j < nmolecule; j++ ) {
        for ( k = 0; mol_list[j][k] >= 0; k++ ) {
            x[mol_list[j][k]] -= dx;
            y[mol_list[j][k]] -= dy;
            z[mol_list[j][k]] -= dz;
        }
    }
}

void wrap_molecules(double *x, double *y, double *z,
                    double a[3], double b[3], double c[3],
                    int mol_list[MAXATOM][MAXATOM], int nmolecule,
                    int bond_list[MAXATOM][MAXATOM], int natom )
/** Wrap the atoms of each molecule to be as close as possible
   to bonded partners. **/
{
    double rx, ry, rz, dx, dy, dz, rx0, ry0, rz0, mtot;
    double mass;
    double x1[2], y1[2], z1[2];
    int wrapped[MAXATOM];
    int j, k, kk, l, ll;

    rx0 = 0.0;
    ry0 = 0.0;
    rz0 = 0.0;
    for ( j = 0; j < natom; j++ ) {
        wrapped[j] = 0;
    }
    for ( j = 0; j < nmolecule; j++ ) {
        rx = 0.0;
        ry = 0.0;
        rz = 0.0;
        mtot = 0.0;
        for ( k = 0; mol_list[j][k] >= 0; k++ ) {
            if ( wrapped[mol_list[j][k]] == 0 ) {
                wrap_pairs(x,y,z,a,b,c,mol_list,nmolecule,bond_list,natom,j,k, wrapped);
            }
        }
    }
}

void wrap_pairs(double *x, double *y, double *z,
                double a[3], double b[3], double c[3],
                int mol_list[MAXATOM][MAXATOM], int nmolecule,
                int bond_list[MAXATOM][MAXATOM], int natom, int j, int k,
                int wrapped[MAXATOM])
/** Wrap atoms bonded to the kth atom of the jth molecule. */
{

    double rx, ry, rz, dx, dy, dz, rx0, ry0, rz0, mtot;
    double mass;
    double x1[2], y1[2], z1[2];
    int kk, l, ll;

    kk = mol_list[j][k];
    for ( l = 0; mol_list[j][l] >= 0; l++ ) {
        ll = mol_list[j][l];
        if ( kk >= natom || ll >= natom ) {
            printf("Error in wrap_molecules: bad index\n");
            exit(1);
        }
        if ( kk != ll && bond_list[kk][ll] == 1 && wrapped[ll] == 0 ) {
            wrap_atom(kk,ll,x,y,z,a,b,c,MAXWRAP,1);
            wrapped[ll] = 1;
            wrap_pairs(x,y,z,a,b,c,mol_list,nmolecule,bond_list,natom,j,l,wrapped);
        }
    }
}

double r2bond(char ea, char eb)
/** Return the square of the bond cutoff distance for
   element pair ea and eb.  This saves having to take the
   square root. **/
{
    int i, j;

    i = ele_index(ea);
    j = ele_index(eb);
    if ( rbond[i][j] < 0.0 ) {
        printf("Error: bond %c-%c was not specified\n", ea, eb);
        exit(1);
    }

    return(rbond[i][j] * rbond[i][j]);
}



int nint(double num)
/* returns the nearest int to the double num */
{

    if (  fabs( num - (int ) num ) <  0.5)
        return( (int ) num );

    else if ( num > 0 )
        return( (int ) num + 1 );
    else
        return( (int ) num - 1 );
}

// NOT USED ANYMORE
//void print_molecule(int mol_list[MAXATOM][MAXATOM],
//                    int bond_list[MAXATOM][MAXATOM],
//                    int *type, char *element, int mol,
//                    int read_charge,
//                    double q[MAXATOM],
//                    double x[MAXATOM], double y[MAXATOM], double z[MAXATOM],
//                    double rcomx[MAXATOM], double rcomy[MAXATOM],
//                    double rcomz[MAXATOM], int printed_atom[MAXATOM],
//                    int dump_structures, int nmolecule, int real_frame, int frame)
//{
//    // Define many variables
//    int k, j, jj, kk, i, a, b, c, unique, retrieved_mol, ok, min_idx, temp, molecule;
//    int bond_count[MAXELE][MAXELE];
//    int ele_count[MAXELE];
//    // Retrieved molecule indices
//    int retrieved_mols[2];
//    // Bond counts for each molecule
//    int bond_counts[MAXATOM][MAXELE][MAXELE];
//    // Element counts for each molecule (used for sorting)
//    int ele_counts[MAXATOM][MAXELE];
//    // Reactant and product molecules
//    int molecules[2][MAXATOM];
//    // Reactant and product molecule counts
//    int molecules_count[2];
//    // List for either mol_list (current or product) or mol_list_prev (past or reactant)
//    int m_list[MAXATOM][MAXATOM];
//    // List for either bond_list (current or product) or bond_list_old (past or reactant
//    int b_list[MAXATOM][MAXATOM];
//    // Array of collisions between molecules and current involved list
//    int collision[MAXATOM][3];
//    // Collision count
//    int collision_count;
//    int atom_involved[MAXATOM];
//    int place;
//    int tester;
//    double charge = 0.0;
//    double dipx, dipy, dipz;
//    double dipole, q0;
//    char ele;
//    char molecule_name[BUFSIZE];
//    char buf[BUFSIZE];
//    int natm_mol;
//    FILE *fmol;
//
//    memset(molecule_name, 0, BUFSIZE);
//    memset(buf, 0, BUFSIZE);
//
//    if ( dump_structures ) {
//            strcat(molecule_name, "molecules/");
//    }
//
//    /** Calculate the overall stoichiometry **/
//    for ( k = 0; k < MAXELE; k++ ) {
//            ele_count[k] = 0;
//            atom_involved[k] = -1;
//    }
//    dipx = 0.0;
//    dipy = 0.0;
//    dipz = 0.0;
//
//    natm_mol = 0;
//    for ( k = 0; mol_list[mol][k] >= 0; k++ ) {
//
//        kk = mol_list[mol][k];
//        q0 = -q[kk];
//        ele = element[type[mol_list[mol][k]]-1];
//        if ( ele == 'C' ) {
//            q0 += 4;
//        } else if ( ele == 'N' ) {
//            q0 += 5;
//        } else if ( ele == 'O' ) {
//            q0 += 6;
//        } else if ( ele == 'F' ) {
//            q0 += 7;
//        } else if ( ele == 'H' ) {
//            q0 += 1;
//        } else {
//            printf("Error: charge for element %c is not known\n", ele);
//        }
//        charge += q0;
//        {
//            double ddx, ddy, ddz;
//            ddx = q0 * (x[kk]-rcomx[mol]);
//            ddy = q0 * (y[kk]-rcomy[mol]);
//            ddz = q0 * (z[kk]-rcomz[mol]);
//
//            dipx += ddx;
//            dipy += ddy;
//            dipz += ddz;
//        }
//				ele_count[type[mol_list[mol][k]]-1]++;
//				natm_mol++;
//		}
//
//		for ( k = 0; mol_list[mol][k] >= 0; k++ ) {
//            kk = mol_list[mol][k];
//            ele = element[type[mol_list[mol][k]]-1];
//            if ( ele == 'C' ) {
//                q0 += 4;
//            } else if ( ele == 'N' ) {
//                q0 += 5;
//            } else if ( ele == 'O' ) {
//                q0 += 6;
//            } else if ( ele == 'F' ) {
//                q0 += 7;
//            } else if ( ele == 'H' ) {
//                q0 += 1;
//            } else {
//                printf("Error: charge for element %c is not known\n", ele);
//            }
//		}
//
//		for ( k = 0; k < MAXELE; k++ ) {
//            if ( ele_count[k] > 0 ) {
//                sprintf(buf, "%c%d ", element[k], ele_count[k]);
//                strncat(molecule_name, buf, BUFSIZE);
//            }
//		}
//
//		/** Calculate the number of bonds **/
//		for ( k = 0; k < MAXELE; k++ ) {
//            for ( j = 0; j < MAXELE; j++ ) {
//                bond_count[j][k] = 0;
//            }
//		}
//		for ( k = 0; mol_list[mol][k] >= 0; k++ ) {
//            kk = mol_list[mol][k];
//            for ( j = 0; mol_list[mol][j] >= 0; j++ ) {
//                jj = mol_list[mol][j];
//                if ( bond_list[kk][jj] == 1 && jj > kk ) {
//                    #if (0)
//                        printf("Bonded: %d-%c %d-%c \n", kk, element[type[kk]-1],
//                       jj, element[type[jj]-1]);
//                    #endif
//                    if ( type[kk] > type[jj] ) {
//                        bond_count[type[kk]-1][type[jj]-1]++;
//                    } else {
//                        bond_count[type[jj]-1][type[kk]-1]++;
//                    }
//                }
//            }
//		}
//
//		for ( k = 0; k < MAXELE; k++ ) {
//            for ( j = 0; j < MAXELE; j++ ) {
//                if ( bond_count[j][k] > 0 ) {
//                    sprintf(buf,"%d(%c-%c) ", bond_count[j][k], element[j], element[k]);
//                    strncat(molecule_name, buf, BUFSIZE);
//                }
//            }
//		}
//    if (count_molecules) {
//        // Adjust lifetimes
//        int found = -1;
//        for (k = 0; k < MAXATOM; k++) {
//            if (strcmp(lifetime_names[real_frame][k], molecule_name) == 0) {
//            found = k;
//            }
//        }
//        if (found >= 0) {
//            lifetime[real_frame][found] += 1;
//        } else {
//            strcpy(lifetime_names[real_frame][lifetime_count[real_frame]], molecule_name);
//            lifetime[real_frame][lifetime_count[real_frame]] = 1;
//            lifetime_count[real_frame]++;
//        }
//    }
//
//    // Reset the molecules array and the counts
//    for (k = 0; k < 2; k++) {
//        for (j = 0; j < MAXATOM; j++) {
//            molecules[k][j] = -1;
//        }
//        molecules_count[k] = 0;
//    }
//
//    // Loop through all bond changes (broken or formed)
//    for (k = 0; bond_changes[mol][k][0] >= 0; k++) {
//        int x, y, z;
//        // Loop through both atoms in each break or formation
//        for (z = 1; z <= 2; z++) {
//            retrieved_mols[0] = retrieve_molecule(bond_changes[mol][k][z], mol_list);
//            retrieved_mols[1] = retrieve_molecule(bond_changes[mol][k][z], mol_list_prev);
//            // This considers reactions with ONE INTERmolecular bond being FORMED
//            if (bonding_count <= 1 && z == 1 && retrieved_mols[1] != retrieve_molecule(bond_changes[mol][k][2], mol_list_prev)) {
//                for (x = 0; x < 2; x++) {
//                    for (a = 0; a < 2; a++) {
//                        local_atoms[retrieved_mols[x]][a] = bond_changes[mol][k][a + 1];
//                    }
//                }
//                bonding_count++;
//            }
//            // Loop through past and present frames (reactants and products)
//            for (x = 0; x < 2; x++) {
//                // Retrieve molecule containing the current atom from corresponding molecule list
//                retrieved_mol = retrieved_mols[x];
//                // If there is a molecule associated with that atom, continue
//                if (retrieved_mol >= 0) {
//                    // Check that it's unique to prevent duplicates in the molecules array
//                    unique = 1;
//                    for (a = 0; molecules[x][a] >= 0; a++) {
//                        if (molecules[x][a] == retrieved_mol) {
//                            unique = 0;
//                        }
//                    }
//                    if (unique) {
//                    // Add that molecule index to the molecules array
//                    molecules[x][molecules_count[x]] = retrieved_mol;
//                    molecules_count[x]++;
//                    }
//                }
//            }
//        }
//    }
//
//    // If there are any changes in molecules detected, continue
//    if (molecules[0][0] >= 0) {
//        // Reset collision array and count
//        for (k = 0; k < MAXATOM; k++) {
//            for (j = 0; j < 3; j++) {
//                collision[k][j] = -1;
//            }
//        }
//        collision_count = 0;
//
//        // Check for collisions with current list of molecules
//        for (k = 0; k < 2; k++) {
//            for (j = 0; molecules[k][j] >= 0; j++) {
//                for (i = 0; involved[i][k][0] >= 0; i++) {
//                    for (a = 0; involved[i][k][a] >= 0; a++) {
//                        // If there is a collision, add it to the collision array
//                        if (involved[i][k][a] == molecules[k][j]) {
//                            collision[collision_count][0] = i;
//                            collision[collision_count][1] = k;
//                            collision[collision_count][2] = a;
//                            collision_count++;
//                        }
//                    }
//                }
//            }
//        }
//        // Check if there was a collision
//        if (collision[0][0] > -1) {
//        // If there was a collision, merge the molecules array into the reaction it collided with
//            for (k = 0; k < 2; k++) {
//                // Do this for past and present
//                for (j = 0; molecules[k][j] >= 0; j++) {
//                    for (i = 0; involved[collision[0][0]][k][i] >= 0; i++) {
//                    // Find the end of the molecule list for the current reaction
//                    }
//                    // Ensure that duplicates are not added
//                    ok = 1;
//                    for (a = 0; collision[a][0] >= 0; a++) {
//                        // Use k instead of collision[a][1] to prevent identifying false duplicates
//                        if (molecules[k][j] == involved[collision[a][0]][k][collision[a][2]]) {
//                            ok = 0;
//                        }
//                    }
//                    // If there are no duplicates, add the molecule to the current reaction
//                    if (ok) {
//                        involved[collision[0][0]][k][i] = molecules[k][j];
//                    }
//                }
//            }
//        } else {
//            // If there was not a collision, add the molecules to a new reaction
//            for (k = 0; k < 2; k++) {
//                for (j = 0; molecules[k][j] >= 0; j++) {
//                    involved[involved_reac_count][k][j] = molecules[k][j];
//                }
//            }
//            involved_reac_count++;
//        }
//    }
//
//    // Print reaction after last molecule in frame with 1 intermolecular bond formed
//    if (mol == nmolecule - 1 && bonding_count == 1 && involved[1][0][0] == -1) { // && involved[1][0][0] == -1
//        // Loop through reactions in involved
//        for (k = 0; involved[k][0][0] >= 0; k++) {
//            // Print buffer
//            char xyz_buffer[1000];
//            int xyz_count = 100; // HACK: issue at frame #254
//            int features[4][9];
//
//            xyz_count += sprintf(xyz_buffer + xyz_count, "\n");
//
//            for (j = 0; j < 4; j++) {
//                for (i = 0; i < 8; i++) {
//                    features[j][i] = -1;
//                }
//            }
//            int features_count = 0;
//
//            int atoms[2] = {-1, -1};
//
//            for (j = 1; j >= 0; j--) {
//                for (i = 0; involved[k][j][i] >= 0; i++) {
//                    int molecule = involved[k][j][i];
//                    if (atoms[0] < 0) {
//                        for (a = 0; a < 2; a++) {
//                            atoms[a] = local_atoms[molecule][a];
//                        }
//                    }
//                }
//            }
//            int total_formed = 0;
//            int total_broken = 0;
//
//            place = 0;
//            for (jj = 0; involved[k][1][jj]>=0; jj++){
//                molecule = involved[k][1][jj];
//                tester = 0;
//                for (a = 0; bond_changes[molecule][a][0] >= 0; a++) {
//                    for (c = 1; c<3; c++){
//                        for (b=0; atom_involved[b]>=0; b++){
//                            if (atom_involved[b] == bond_changes[molecule][a][c]){
//                                tester = 1;
//                            }
//                        }
//                        if (tester == 0){
//                            atom_involved[place] = bond_changes[molecule][a][c];
//                            place++;
//                        }
//                        tester = 0;
//                    }
//                }
//            }
//
//            int featfeat[2*place][5];
//            int ll;
//            for (jj=0;jj<2*place;jj++){
//                for (kk=0;kk<5; kk++){
//                    featfeat[jj][kk] = 50;
//                }
//            }
//            for (ll = 0; ll<place; ll++){
//                featfeat[ll][0] = type[atom_involved[ll]];
//                featfeat[ll+place][0] = type[atom_involved[ll]];
//                featfeat[ll][4] = atom_involved[ll];
//                featfeat[ll+place][4] = atom_involved[ll];
//            }
//
//
//            for (j=0; j<2; j++){
//                if (j == 0) {
//                    memcpy(m_list, mol_list, sizeof m_list);
//                    memcpy(b_list, bond_list, sizeof b_list);
//                    features_count = 2;
//                } else {
//                    memcpy(m_list, mol_list_prev, sizeof m_list);
//                    memcpy(b_list, bond_prev, sizeof b_list);
//                }
//                for (ll=0; ll<place; ll++){
//                    int atom = atom_involved[ll];
//                    molecule = -1;
//                    for (a=0;involved[k][j][a]>=0; a++){
//                        if(is_in_molecule(atom, involved[k][j][a], m_list)){
//                            molecule = involved[k][j][a];
//                        }
//                    }
//                    if(molecule == -1){
//                        printf("Error: an atom is not in involved molecule");
//                    }
//
//                    for ( b = 0; b < MAXELE; b++ ) {
//                        for ( c = 0; c < MAXELE; c++ ) {
//                            bond_count[b][c] = 0;
//                        }
//                    }
//
//                    for ( b = 0; m_list[molecule][b] >= 0; b++ ) {
//                        jj = m_list[molecule][b];
//                        if ( b_list[atom][jj] == 1) {
//                            if ( type[atom] > type[jj] ) {
//                                bond_count[type[atom]-1][type[jj]-1]++;
//                            } else {
//                                bond_count[type[jj]-1][type[atom]-1]++;
//                            }
//                        }
//                    }
//
//                    int local_counts[3] = {0, 0, 0}; // (C-C), (H-H), (H-C)
//
//                    for ( b = 0; b < MAXELE; b++ ) {
//                        for ( c = 0; c < MAXELE; c++ ) {
//                            if ( bond_count[c][b] > 0 ) {
//                                if (c == 0 && b == 0) local_counts[0] += bond_count[c][b];
//                                if (c == 1 && b == 1) local_counts[1] += bond_count[c][b];
//                                if (c == 1 && b == 0) local_counts[2] += bond_count[c][b];
//                            }
//                        }
//                    }
//                    for(kk=0;kk<3;kk++){
//                        featfeat[ll+(1-j)*place][kk+1]= local_counts[kk];
//                    }
//                }
//            }
//            int sort_featfeat[2*place][5];
//            int sort_help[place][3];
//            int indexes[place];
//            for(jj=0; jj<place; jj++){
//                sort_help[jj][0] = 1*featfeat[jj][3] + 10*featfeat[jj][2] + 100*featfeat[jj][1] + 1000*featfeat[jj][0];
//                sort_help[jj][1] = jj;
//                sort_help[jj][2] = 1*featfeat[jj+place][3] + 10*featfeat[jj+place][2] + 100*featfeat[jj+place][1] + 1000*featfeat[jj+place][0];
//            }
//
//            for(jj=0; jj<place; jj++){
//		        for(kk=0; kk<place-1; kk++){
//			        if(sort_help[kk][0]>sort_help[kk+1][0]){
//                        int temp = sort_help[kk+1][0];
//                        sort_help[kk+1][0] = sort_help[kk][0];
//                        sort_help[kk][0] = temp;
//                        temp = sort_help[kk+1][1];
//                        sort_help[kk+1][1] = sort_help[kk][1];
//                        sort_help[kk][1] = temp;
//                        temp = sort_help[kk+1][2];
//                        sort_help[kk+1][2] = sort_help[kk][2];
//                        sort_help[kk][2] = temp;
//			        }
//			        else if(sort_help[kk][0]==sort_help[kk+1][0]) {
//			            if (sort_help[kk][2]>sort_help[kk+1][2]){
//                            int temp = sort_help[kk+1][0];
//                            sort_help[kk+1][0] = sort_help[kk][0];
//                            sort_help[kk][0] = temp;
//                            temp = sort_help[kk+1][1];
//                            sort_help[kk+1][1] = sort_help[kk][1];
//                            sort_help[kk][1] = temp;
//                            temp = sort_help[kk+1][2];
//                            sort_help[kk+1][2] = sort_help[kk][2];
//                            sort_help[kk][2] = temp;
//			            }
//			        }
//		        }
//	        }
//
//
//	        for(jj=0;jj<place;jj++){
//	            for(kk=0;kk<5;kk++){
//	                sort_featfeat[jj][kk] = featfeat[sort_help[jj][1]][kk];
//	                sort_featfeat[jj+place][kk] = featfeat[sort_help[jj][1]+place][kk];
//	            }
//	        }
//
//	        for(jj=0;jj<2*place;jj++){
//	            for(kk=0;kk<5;kk++){
//	                featfeat[jj][kk] = sort_featfeat[jj][kk];
//	            }
//	        }
//
//            printf("Reaction with atoms: ");
//            for (ll=0; ll<2*place;ll++){
//                for(jj=0; jj<4; jj++){
//                    printf("%d ", featfeat[ll][jj]);
//                }
//                if (ll == place - 1){
//                    printf("=> ");
//                }
//                else if (ll == 2*place - 1){
//                    printf("\n");
//                }
//                else{
//                    printf("+ ");
//                }
//            }
//            printf("Atoms involved: ");
//            for (ll=0; ll<place;ll++){
//                printf("%d ", featfeat[ll][4]);
//            }
//            printf("\n");
//
//            printf("Reaction with molecules: ");
//            for (j = 1; j >= 0; j--) {
//                // Select correct molecule and bond lists based on reactant or product
//                if (j == 0) {
//                    memcpy(m_list, mol_list, sizeof m_list);
//                    memcpy(b_list, bond_list, sizeof b_list);
//                    features_count = 2;
//                } else {
//                    memcpy(m_list, mol_list_prev, sizeof m_list);
//                    memcpy(b_list, bond_prev, sizeof b_list);
//                }
//
//                // Get element counts and bond counts for each molecule in reactants/products
//                for (i = 0; involved[k][j][i] >= 0; i++) {
//                    int molecule = involved[k][j][i];
//
//                    for ( a = 0; a < MAXELE; a++ ) {
//                        ele_counts[i][a] = 0;
//        		    }
//
//                    for ( a = 0; m_list[molecule][a] >= 0; a++ ) {
//                        ele_counts[i][type[m_list[molecule][a]]-1]++;
//                    }
//
//                    // Reset bond counts
//                    for ( a = 0; a < MAXELE; a++ ) {
//        				for ( b = 0; b < MAXELE; b++ ) {
//                            bond_counts[i][a][b] = 0;
//        				}
//        	    	}
//
//                    // Calculate bond counts
//                    for ( a = 0; m_list[molecule][a] >= 0; a++ ) {
//        				kk = m_list[molecule][a];
//        				for ( b = 0; m_list[molecule][b] >= 0; b++ ) {
//                            jj = m_list[molecule][b];
//                            if ( b_list[kk][jj] == 1 && jj > kk ) {
//                                if ( type[kk] > type[jj] ) {
//                                    bond_counts[i][type[kk]-1][type[jj]-1]++;
//                                } else {
//                                    bond_counts[i][type[jj]-1][type[kk]-1]++;
//                                }
//                            }
//        				}
//            		}
//
//                    if (j == 1) {
//                        for (a = 0; bond_changes[molecule][a][0] >= 0; a++) {
//                            if (bond_changes[molecule][a][0] == 0) {
//                                total_broken += 1;
//                            } else {
//                                total_formed += 1;
//                            }
//                        }
//                    }
//                }
//
//                // Sort molecules with a selection sort in reactants/products based on element count to prevent different permutations of the same reaction
//                for (i = 0; involved[k][j][i + 1] >= 0; i++) {
//                    min_idx = i;
//                    for (a = i + 1; involved[k][j][a] >= 0; a++) {
//                        if ((ele_counts[a][1] < ele_counts[min_idx][1] && ele_counts[a][0] > 0) || (ele_counts[a][0] > 0 && ele_counts[min_idx][0] == 0)) {
//                            min_idx = a;
//                        }
//                    }
//
//                    // Adjust not only the element counts, but also the array of molecules and the bond counts
//                    temp = involved[k][j][min_idx];
//                    involved[k][j][min_idx] = involved[k][j][i];
//                    involved[k][j][i] = temp;
//                    for (a = 0; a < MAXELE; a++) {
//                        temp = ele_counts[min_idx][a];
//                        ele_counts[min_idx][a] = ele_counts[i][a];
//                        ele_counts[i][a] = temp;
//                        for (b = 0; b < MAXELE; b++) {
//                            temp = bond_counts[min_idx][a][b];
//                            bond_counts[min_idx][a][b] = bond_counts[i][a][b];
//                            bond_counts[i][a][b] = temp;
//                        }
//                    }
//                }
//
//                // BEGINNING NEW CODE
//
//                // Loop through each molecule to print it
//                for (i = 0; involved[k][j][i] >= 0; i++) {
//                    // Retrieve the molecule index
//                    molecule = involved[k][j][i];
//
//                    if (j == 0 || involved[k][j][i + 1] == -1) {
//                        xyz_count += sprintf(xyz_buffer + xyz_count, "\n");
//                    }
//
//                    // Print out the elements and their counts
//        		    for ( a = 0; a < MAXELE; a++ ) {
//        				if ( ele_counts[i][a] > 0 ) {
//                            printf("%c%d ", element[a], ele_counts[i][a]);
//                            xyz_count += sprintf(xyz_buffer + xyz_count, "%c%d ", element[a], ele_counts[i][a]);
//        				}
//        		    }
//
//                    features[features_count][0] = ele_counts[i][0]; // carbon count
//                    features[features_count][1] = ele_counts[i][1]; // hydrogen count
//
//                    int counts[3] = {0, 0, 0}; // (C-C), (H-H), (H-C)
//
//                    // Print bond counts and types
//                    for ( a = 0; a < MAXELE; a++ ) {
//        				for ( b = 0; b < MAXELE; b++ ) {
//                            if ( bond_counts[i][b][a] > 0 ) {
//                                printf("%d(%c-%c) ", bond_counts[i][b][a], element[b], element[a]);
//                                xyz_count += sprintf(xyz_buffer + xyz_count, "%d(%c-%c) ", bond_counts[i][b][a], element[b], element[a]);
//                                if (b == 0 && a == 0) counts[0] += bond_counts[i][b][a];
//                                if (b == 1 && a == 1) counts[1] += bond_counts[i][b][a];
//                                if (b == 1 && a == 0) counts[2] += bond_counts[i][b][a];
//                            }
//        				}
//        		    }
//
//                    xyz_count += sprintf(xyz_buffer + xyz_count, "\n");
//
//                    features[features_count][2] = counts[0]; // C-C bond count
//                    features[features_count][3] = counts[1]; // H-H bond count
//
//                    int formal_charge = 0;
//
//                    // Calculate formal charge
//                    for ( a = 0; m_list[molecule][a] >= 0; a++ ) {
//                        int formal_counts[3] = {0, 0, 0}; // (C-C), (H-H), (H-C)
//                        for ( b = 0; b < MAXELE; b++ ) {
//                            for ( c = 0; c < MAXELE; c++ ) {
//                                bond_count[b][c] = 0;
//                            }
//                        }
//                        kk = m_list[molecule][a];
//                        for ( b = 0; m_list[molecule][b] >= 0; b++ ) {
//                            jj = m_list[molecule][b];
//                            if ( b_list[kk][jj] == 1) {
//                                if ( type[kk] > type[jj] ) {
//                                    bond_count[type[kk]-1][type[jj]-1]++;
//                                } else {
//                                    bond_count[type[jj]-1][type[kk]-1]++;
//                                }
//                            }
//                        }
//                        for ( b = 0; b < MAXELE; b++ ) {
//                            for ( c = 0; c < MAXELE; c++ ) {
//                                if ( bond_count[c][b] > 0 ) {
//                                    if (c == 0 && b == 0) formal_counts[0] += bond_count[c][b];
//                                    if (c == 1 && b == 1) formal_counts[1] += bond_count[c][b];
//                                    if (c == 1 && b == 0) formal_counts[2] += bond_count[c][b];
//                                }
//                            }
//                        }
//                        if (type[kk] - 1 == 0) {
//                            formal_charge += 4 - (2 * (4 - (formal_counts[0] + formal_counts[2]))) - (formal_counts[0] + formal_counts[2]);
//                        } else {
//                            formal_charge += 1 - (2 * (1 - (formal_counts[1] + formal_counts[2]))) - (formal_counts[1] + formal_counts[2]);
//                        }
//                    }
//
//                    features[features_count][4] = formal_charge; // sum of formal charges
//
//                    int is_contained[2];
//                    for (a = 0; a < 2; a++) {
//                        is_contained[a] = is_in_molecule(atoms[a], molecule, m_list);
//                    }
//                    int loop = is_contained[0] && is_contained[1] ? 2 : (!is_contained[0] && !is_contained[1] ? 0 : 1);
//                    for (a = 0; a < loop; a++) {
//                        int atom = loop == 2 ? atoms[a] : (is_contained[0] ? atoms[0] : atoms[1]);
//                        if (a == 1) {
//                            for (b = 0; b < 5; b++) {
//                                features[features_count][b] = features[features_count - 1][b];
//                            }
//                        }
//
//                        for ( b = 0; b < MAXELE; b++ ) {
//                            for ( c = 0; c < MAXELE; c++ ) {
//                                bond_count[b][c] = 0;
//                            }
//                        }
//
//                        for ( b = 0; m_list[molecule][b] >= 0; b++ ) {
//                            jj = m_list[molecule][b];
//                            if ( b_list[atom][jj] == 1) {
//                                if ( type[atom] > type[jj] ) {
//                                    bond_count[type[atom]-1][type[jj]-1]++;
//                                } else {
//                                    bond_count[type[jj]-1][type[atom]-1]++;
//                                }
//                            }
//                        }
//
//                        int local_counts[3] = {0, 0, 0}; // (C-C), (H-H), (H-C)
//
//                        for ( b = 0; b < MAXELE; b++ ) {
//                            for ( c = 0; c < MAXELE; c++ ) {
//                                if ( bond_count[c][b] > 0 ) {
//                                    if (c == 0 && b == 0) local_counts[0] += bond_count[c][b];
//                                    if (c == 1 && b == 1) local_counts[1] += bond_count[c][b];
//                                    if (c == 1 && b == 0) local_counts[2] += bond_count[c][b];
//                                }
//                            }
//                        }
//
//                        int atom_type = type[atom] - 1;
//
//                        // Octet rule calculation
//                        if (atom_type == 0) {
//                            features[features_count][5] = (local_counts[0] + local_counts[2]) - 4;
//                        } else {
//                            features[features_count][5] = (local_counts[1] + local_counts[2]) - 1;
//                        }
//
//                        features[features_count][6] = local_counts[0]; // C-C bond count
//                        features[features_count][7] = local_counts[2]; // H-C bond count
//
//                        features[features_count][8] = atom_type; // atom type
//
//                        features_count++;
//                    }
//
//                    // Print plus only if there is another reactant/product
//                    if (involved[k][j][i + 1] >= 0) {
//                        printf("+ ");
//                    }
//                }
//
//                // Print arrow only if products are next
//                if (j == 1) {
//                    printf("=> ");
//                    xyz_count += sprintf(xyz_buffer + xyz_count, "\n----\n");
//                }
//            }
//            printf("\n");
//        }
//    }
//
//    if ( dump_structures ) {
//        strncat(molecule_name, ".xyz", BUFSIZE);
//        fmol = fopen(molecule_name, "w");
//        if ( fmol == NULL ) {
//            printf("Error: could not open file %s\n", molecule_name);
//            exit(1);
//        }
//        fprintf(fmol,"%d\n\n", k);
//        for ( k = 0; mol_list[mol][k] >= 0; k++ ) {
//            j = mol_list[mol][k];
//            ele = element[type[mol_list[mol][k]]-1];
//            fprintf(fmol,"%c %11.4f %11.4f %11.4f\n", ele, x[j], y[j], z[j]);
//        }
//        fclose(fmol);
//    }
//
//    /** 4.8 is a unit conversion factor. */
//    dipole = 4.8 * sqrt(dipx * dipx + dipy * dipy + dipz * dipz);
//    if ( read_charge ) {
//        printf("   Charge = %f e\n", charge);
//        printf("   Dipole moment = %f D\n", dipole);
//    }
//
//}


static void wrap_in_box(double a[3], double b[3], double c[3],
                        double invbox[3][3],
                        int natom, double x[MAXATOM], double y[MAXATOM],
                        double z[MAXATOM])
/** This function wraps all the atoms into the primitive simulation
   box. **/
{
    int j;
    double ca, cb, cc;

    for ( j = 0; j < natom; j++ ) {
        ca = invbox[0][0] * x[j] + invbox[1][0] * y[j] + invbox[2][0] * z[j];
        cb = invbox[0][1] * x[j] + invbox[1][1] * y[j] + invbox[2][1] * z[j];
        cc = invbox[0][2] * x[j] + invbox[1][2] * y[j] + invbox[2][2] * z[j];

        ca -= nint(ca);
        cb -= nint(cb);
        cc -= nint(cc);

        x[j] = ca * a[0] + cb * b[0] + cc * c[0];
        y[j] = ca * a[1] + cb * b[1] + cc * c[1];
        z[j] = ca * a[2] + cb * b[2] + cc * c[2];

    }
}

void inversebox(double a[3], double b[3], double c[3],
                double invbox[3][3])
/** Calculate the inverse box maxtrix, used for finding
   crystal coordinates. */
{
    double xhlp;
    xhlp=-a[2]*b[1]*c[0]+ a[1]*b[2]*c[0];
    xhlp=xhlp+a[2]*b[0]*c[1];
    xhlp=xhlp-a[0]*b[2]*c[1];
    xhlp=xhlp-a[1]*b[0]*c[2];
    xhlp=xhlp+a[0]*b[1]*c[2];

    invbox[0][0]=(-b[2]*c[1]+b[1]*c[2])/xhlp;
    invbox[1][0]=(b[2]*c[0]-b[0]*c[2])/xhlp;
    invbox[2][0]=(-b[1]*c[0]+b[0]*c[1])/xhlp;

    invbox[0][1]=(a[2]*c[1]-a[1]*c[2])/xhlp;
    invbox[1][1]=(-a[2]*c[0]+a[0]*c[2])/xhlp;
    invbox[2][1]=(a[1]*c[0]-a[0]*c[1])/xhlp;

    invbox[0][2]=(-a[2]*b[1]+a[1]*b[2])/xhlp;
    invbox[1][2]=(a[2]*b[0]-a[0]*b[2])/xhlp;
    invbox[2][2]=(-a[1]*b[0]+a[0]*b[1])/xhlp;
}


double atm_mass(char ea)
/* Return the atomic mass of the given element. */
{
    if ( ea == 'H' ) {
        return(1.0);
    } else if ( ea == 'C' ) {
        return(12.0);
    } else if ( ea == 'N' ) {
        return(14.0);
    } else if ( ea == 'O' ) {
        return(16.0);
    } else if ( ea == 'F' ) {
        return(19.0);
    } else {
        printf("Error: element %d is unknown\n", ea );
        exit(1);
    }
    return(1.0);
}

void read_bonds(double rbond[MAXELE][MAXELE], int *duration, int *xyz_copies, double *time_step,
                int *dump_structures)
/** Read in bond distances and other parameters.
   rbond - The maximum bond distance.
   duration - The bond lifetime in simulation steps.
   xyz_copies - The number of copies in each direction for the xyz output file.
   time_step - The time between each frame stored.
   dump_structures - Whether to write out individual molecule structure files.
 **/
{
    FILE *fbond;
    char ele1, ele2;
    char buf1[BUFSIZ], buf2[BUFSIZ];
    double r;
    int i, j, k;

    for ( i = 0; i < MAXELE; i++ ) {
        for ( j = 0; j < MAXELE; j++ ) {
            /** Put in a negative value so that unspecified bonds can be caught. */
            rbond[i][j] = -1.0;
        }
    }

    fbond = fopen("bonds.dat","r");
    if ( fbond == NULL ) {
        printf("Error: could not open bonds.dat\n");
        exit(1);
    }
    fscanf(fbond, "%d", duration);
    printf("Bonds require a lifetime of %d steps\n", *duration);

    fscanf(fbond, "%d", xyz_copies);
    printf("%d replicas of the simulation cell will be propagated in each direction\n",
       *xyz_copies);

    fscanf(fbond, "%lf\n", time_step);
    printf("The time between frames read in = %11.7e seconds\n", *time_step);
    if ( *time_step < 0.0 || *time_step > 1.0e-10 ) {
        printf("The time step was given as %21.14e seconds.  This seems strange\n", *time_step);
        exit(1);
    }

    fscanf(fbond, "%d", dump_structures);
    if ( *dump_structures != 0 ) {
        printf("Individual molecule xyz files will be created in the molecules directory.\n");
    } else {
        printf("Individual molecule xyz file will not be created.\n");
    }
    printf("Bond distances:\n");

    do {
        if ( fscanf(fbond, "%s %s %lf", buf1, buf2, &r) != 3 ) {
            break;
        }
        ele1 = buf1[0];
        ele2 = buf2[0];
        i = ele_index(ele1);
        j = ele_index(ele2);
        if ( j < i ) {
            k = i;
            i = j;
            j = k;
        }
        printf("%c %c %f\n", ele1, ele2, r);
        if ( rbond[i][j] > 0.0 || rbond[j][i] > 0.0 ) {
            printf("Error: Bond distance for %c %c has already been specified\n",
               ele1, ele2);
            exit(1);
        }
        rbond[i][j] = r;

    } while ( !feof(fbond) && !ferror(fbond));
    fclose(fbond);

    for ( i = 0; i < MAXELE; i++ ) {
        for ( j = i + 1; j < MAXELE; j++ ) {
            rbond[j][i] = rbond[i][j];
        }
    }
}

int ele_index(char ea)
/** Returns the index of the given element */
{
    int i;
    if ( ea == 'H' ) {
        i = 0;
    } else if ( ea == 'C' ) {
        i = 1;
    } else if ( ea == 'N' ) {
        i = 2;
    } else if ( ea == 'O' ) {
        i = 3;
    } else if ( ea == 'F' ) {
        i = 4;
    } else {
        printf("Error: element %d is unknown\n", ea );
        exit(1);
    }
    return i;
}

double boxvol(double a[3], double b[3], double c[3])
/* Calculates the volume of the simulation box, given the
   box cell vectors a, b, and c. */
{

    double bv;
    double xhlp;
    double ax, ay, az, bx, by, bz, cx, cy, cz;

    ax = a[0];
    ay = a[1];
    az = a[2];
    bx = b[0];
    by = b[1];
    bz = b[2];
    cx = c[0];
    cy = c[1];
    cz = c[2];

    xhlp=-az*by*cx+ ay*bz*cx;
    xhlp=xhlp+az*bx*cy;
    xhlp=xhlp-ax*bz*cy;
    xhlp=xhlp-ay*bx*cz;
    xhlp=xhlp+ax*by*cz;

    bv = fabs(xhlp);
    return bv;
}

int retrieve_molecule(int atom_id, int mol_list[MAXATOM][MAXATOM])
{
    int x, y;
    for (x = 0; mol_list[x][0] >= 0; x++) {
        for (y = 0; mol_list[x][y] >= 0; y++) {
            if (mol_list[x][y] == atom_id) {
                return x;
            }
        }
    }
    return -1;
}

int is_in_molecule(int atom_id, int mol_id, int mol_list[MAXATOM][MAXATOM])
{
    int x;
    for (x = 0; mol_list[mol_id][x] >= 0; x++) {
        if (mol_list[mol_id][x] == atom_id) {
            return 1;
        }
    }
    return 0;
}

