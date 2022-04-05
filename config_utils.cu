#include "globals.h"
void allocate_host_particles(void);

// Must know Dim before read_config() is called! //
void read_config(const char *nm) {

	int i, di, ind, ltp;
	double dx, dy;
	char tt[120];

	FILE* inp;
	inp = fopen(nm, "r");
	if (inp == NULL) {
		char death[50];
		sprintf(death, "Failed to open %s!\n", nm);
		die(death);
	}

	(void)!fgets(tt, 120, inp);
	(void)!fgets(tt, 120, inp);

	(void)!fscanf(inp, "%d", &ns);      (void)!fgets(tt, 120, inp);
	(void)!fscanf(inp, "%d", &n_total_bonds);  (void)!fgets(tt, 120, inp);
	(void)!fscanf(inp, "%d", &n_total_angles);  (void)!fgets(tt, 120, inp);

	(void)!fgets(tt, 120, inp);

	(void)!fscanf(inp, "%d", &ntypes);  (void)!fgets(tt, 120, inp);
	(void)!fscanf(inp, "%d", &nbond_types);  (void)!fgets(tt, 120, inp);
	(void)!fscanf(inp, "%d", &nangle_types);  (void)!fgets(tt, 120, inp);

	// Read in box shape
	(void)!fgets(tt, 120, inp);
	(void)!fscanf(inp, "%lf %lf", &dx, &dy);   (void)!fgets(tt, 120, inp);
	L[0] = (float)(dy - dx);
	(void)!fscanf(inp, "%lf %lf", &dx, &dy);   (void)!fgets(tt, 120, inp);
	L[1] = (float)(dy - dx);
	(void)!fscanf(inp, "%lf %lf", &dx, &dy);   (void)!fgets(tt, 120, inp);
	if (Dim > 2)
		L[2] = (float)(dy - dx);
	else
		L[2] = 1.0f;

	V = 1.0f;
	for (i = 0; i < Dim; i++) {
		Lh[i] = 0.5f * L[i];
		V *= L[i];
	}


	// Allocate memory for positions //
	allocate_host_particles();
	
	printf("Particle memory allocated on host!\n");

	(void)!fgets(tt, 120, inp);

	// Read in particle masses
	(void)!fgets(tt, 120, inp);
	(void)!fgets(tt, 120, inp);
	for (i = 0; i < ntypes; i++) {
		(void)!fscanf(inp, "%d %lf", &di, &dx); (void)!fgets(tt, 120, inp);
		mass[di - 1] = (float) dx;
	}
	(void)!fgets(tt, 120, inp);

	
	// Read in atomic positions
	(void)!fgets(tt, 120, inp);
	(void)!fgets(tt, 120, inp);


	for (i = 0; i < ns; i++) {
		if (feof(inp)) die("Premature end of input.conf!");

		(void)!fscanf(inp, "%d %d %d", &ind, &di, &ltp);
		ind -= 1;

		molecID[ind] = di - 1;
		tp[ind] = ltp - 1;

		for (int j = 0; j < Dim; j++) {
			(void)!fscanf(inp, "%lf", &dx);
			x[ind][j] = (float) dx;
		}

		(void)!fgets(tt, 120, inp);
	}
	(void)!fgets(tt, 120, inp);


	// Read in bond information
	(void)!fgets(tt, 120, inp);
	(void)!fgets(tt, 120, inp);
	
	for (i = 0; i < ns; i++)
		n_bonds[i] = 0;

	for (i = 0; i < n_total_bonds; i++) {
		if (feof(inp)) die("Premature end of input.conf!");
		(void)!fscanf(inp, "%d", &di);
		(void)!fscanf(inp, "%d", &di);
		int b_type = di - 1;

		(void)!fscanf(inp, "%d", &di);
		int i1 = di - 1;

		(void)!fscanf(inp, "%d", &di);
		int i2 = di - 1;

		if (i2 < i1) {
			di = i2;
			i2 = i1;
			i1 = di;
		}

		bonded_to[i1][n_bonds[i1]] = i2;
		bond_type[i1][n_bonds[i1]] = b_type;
		n_bonds[i1]++;

		bonded_to[i2][n_bonds[i2]] = i1;
		bond_type[i2][n_bonds[i2]] = b_type;
		n_bonds[i2]++;

	}
	(void)!fgets(tt, 120, inp);



	// Read in angle information
	for (i = 0; i < ns; i++)
		n_angles[i] = 0;

	(void)!fgets(tt, 120, inp);
	(void)!fgets(tt, 120, inp);
	for (i = 0; i < n_total_angles; i++) {

		(void)!fscanf(inp, "%d", &di);
		(void)!fscanf(inp, "%d", &di);

		int a_type = di - 1;

		(void)!fscanf(inp, "%d", &di);
		int i1 = di - 1;

		(void)!fscanf(inp, "%d", &di);
		int i2 = di - 1;

		(void)!fscanf(inp, "%d", &di);
		int i3 = di - 1;

		if (i3 < i1) {
			di = i3;
			i3 = i1;
			i1 = di;
		}

		int na = n_angles[i1];
		angle_first[i1][na] = i1;
		angle_mid[i1][na] = i2;
		angle_end[i1][na] = i3;
		angle_type[i1][na] = a_type;
		n_angles[i1] += 1;

		na = n_angles[i2];
		angle_first[i2][na] = i1;
		angle_mid[i2][na] = i2;
		angle_end[i2][na] = i3;
		angle_type[i2][na] = a_type;
		n_angles[i2] += 1;

		na = n_angles[i3];
		angle_first[i3][na] = i1;
		angle_mid[i3][na] = i2;
		angle_end[i3][na] = i3;
		angle_type[i3][na] = a_type;
		n_angles[i3] += 1;

		(void)!fgets(tt, 120, inp);
		
	}
	fclose(inp);
}

void read_charge_config(const char* nm) {

	int i, di, ind, ltp;
	float dx, dy, dcharge;
	char tt[120];

	FILE* inp;
	inp = fopen(nm, "r");
	if (inp == NULL) {
		char death[50];
		sprintf(death, "Failed to open %s!\n", nm);
		die(death);
	}

	(void)!fgets(tt, 120, inp);
	(void)!fgets(tt, 120, inp);

	(void)!fscanf(inp, "%d", &ns);      (void)!fgets(tt, 120, inp);
	(void)!fscanf(inp, "%d", &n_total_bonds);  (void)!fgets(tt, 120, inp);
	(void)!fscanf(inp, "%d", &n_total_angles);  (void)!fgets(tt, 120, inp);

	(void)!fgets(tt, 120, inp);

	(void)!fscanf(inp, "%d", &ntypes);  (void)!fgets(tt, 120, inp);
	(void)!fscanf(inp, "%d", &nbond_types);  (void)!fgets(tt, 120, inp);
	(void)!fscanf(inp, "%d", &nangle_types);  (void)!fgets(tt, 120, inp);

	// Read in box shape
	(void)!fgets(tt, 120, inp);
	(void)!fscanf(inp, "%f %f", &dx, &dy);   (void)!fgets(tt, 120, inp);
	L[0] = (float)(dy - dx);
	(void)!fscanf(inp, "%f %f", &dx, &dy);   (void)!fgets(tt, 120, inp);
	L[1] = (float)(dy - dx);
	(void)!fscanf(inp, "%f %f", &dx, &dy);   (void)!fgets(tt, 120, inp);
	if (Dim > 2)
		L[2] = (float)(dy - dx);
	else
		L[2] = 1.0f;

	V = 1.0f;
	for (i = 0; i < Dim; i++) {
		Lh[i] = 0.5f * L[i];
		V *= L[i];
	}


	// Allocate memory for positions //
	allocate_host_particles();

	printf("Particle memory allocated on host!\n");

	(void)!fgets(tt, 120, inp);

	// Read in particle masses
	(void)!fgets(tt, 120, inp);
	(void)!fgets(tt, 120, inp);
	for (i = 0; i < ntypes; i++) {
		(void)!fscanf(inp, "%d %f", &di, &dx); (void)!fgets(tt, 120, inp);
		mass[di - 1] = float(dx);
	}
	(void)!fgets(tt, 120, inp);


	// Read in atomic positions
	(void)!fgets(tt, 120, inp);
	(void)!fgets(tt, 120, inp);


	for (i = 0; i < ns; i++) {
		if (feof(inp)) die("Premature end of input.conf!");

		(void)!fscanf(inp, "%d %d %d", &ind, &di, &ltp);
		ind -= 1;

		molecID[ind] = di - 1;
		tp[ind] = ltp - 1;

		(void)!fscanf(inp, "%f", &dcharge);
		charges[ind] = dcharge;

		for (int j = 0; j < Dim; j++) {
			(void)!fscanf(inp, "%f", &dx);
			x[ind][j] = dx;
		}

		(void)!fgets(tt, 120, inp);
	}
	(void)!fgets(tt, 120, inp);

	// Read in bond information
	(void)!fgets(tt, 120, inp);
	(void)!fgets(tt, 120, inp);

	for (i = 0; i < ns; i++)
		n_bonds[i] = 0;

	for (i = 0; i < n_total_bonds; i++) {
		(void)!fscanf(inp, "%d", &di);
		(void)!fscanf(inp, "%d", &di);
		int b_type = di - 1;

		(void)!fscanf(inp, "%d", &di);
		int i1 = di - 1;

		(void)!fscanf(inp, "%d", &di);
		int i2 = di - 1;

		if (i2 < i1) {
			di = i2;
			i2 = i1;
			i1 = di;
		}

		bonded_to[i1][n_bonds[i1]] = i2;
		bond_type[i1][n_bonds[i1]] = b_type;
		n_bonds[i1]++;

		bonded_to[i2][n_bonds[i2]] = i1;
		bond_type[i2][n_bonds[i2]] = b_type;
		n_bonds[i2]++;

	}
	(void)!fgets(tt, 120, inp);



	// Read in angle information
	for (i = 0; i < ns; i++)
		n_angles[i] = 0;

	(void)!fgets(tt, 120, inp);
	(void)!fgets(tt, 120, inp);
	for (i = 0; i < n_total_angles; i++) {

		(void)!fscanf(inp, "%d", &di);
		(void)!fscanf(inp, "%d", &di);

		int a_type = di - 1;

		(void)!fscanf(inp, "%d", &di);
		int i1 = di - 1;

		(void)!fscanf(inp, "%d", &di);
		int i2 = di - 1;

		(void)!fscanf(inp, "%d", &di);
		int i3 = di - 1;

		if (i3 < i1) {
			di = i3;
			i3 = i1;
			i1 = di;
		}

		int na = n_angles[i1];
		angle_first[i1][na] = i1;
		angle_mid[i1][na] = i2;
		angle_end[i1][na] = i3;
		angle_type[i1][na] = a_type;
		n_angles[i1] += 1;

		na = n_angles[i2];
		angle_first[i2][na] = i1;
		angle_mid[i2][na] = i2;
		angle_end[i2][na] = i3;
		angle_type[i2][na] = a_type;
		n_angles[i2] += 1;

		na = n_angles[i3];
		angle_first[i3][na] = i1;
		angle_mid[i3][na] = i2;
		angle_end[i3][na] = i3;
		angle_type[i3][na] = a_type;
		n_angles[i3] += 1;

		(void)!fgets(tt, 120, inp);

	}
	fclose(inp);
}

void write_lammps_traj() {

	FILE* otp;
	int i, j;
	if (step == 0)
		otp = fopen(dump_name.c_str(), "w");
	else
		otp = fopen(dump_name.c_str(), "a");

	fprintf(otp, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\n", step, ns);
	fprintf(otp, "ITEM: BOX BOUNDS pp pp pp\n");
	fprintf(otp, "%f %f\n%f %f\n%f %f\n", 0.f, L[0],
		0.f, L[1],
		(Dim == 3 ? 0.f : 1.f), (Dim == 3 ? L[2] : 1.f));

	if ( do_charges )
		fprintf(otp, "ITEM: ATOMS id type mol x y z q\n");
	else
		fprintf(otp, "ITEM: ATOMS id type mol x y z\n");

	for (i = 0; i < ns; i++) {
		fprintf(otp, "%d %d %d  ", i + 1, tp[i] + 1, molecID[i] + 1);
		for (j = 0; j < Dim; j++)
			fprintf(otp, "%f ", x[i][j]);

		for (j = Dim; j < 3; j++)
			fprintf(otp, "%f", 0.f);

		if ( do_charges )
			fprintf(otp, " %f", charges[i]);

		fprintf(otp, "\n");
	}
	fclose(otp);
}

/*
void write_conf_file() {
    FILE* otp;
    int i, j;
    otp = fopen(final_conf_name.c_str(), "w");

    fprintf(otp, "GPU TILD\n\n");
    
    fprintf(otp, "%d atoms\n", ns);
    fprintf(otp, "%d bonds\n", n_total_bonds);
    fprintf(otp, "%d angles\n", n_total_angles);
    fprintf(otp, "\n");

    fprintf(otp, "%d atom types\n", ntypes);
    fprintf(otp, "%d bond types\n", nbond_types);
    fprintf(otp, "%d angle types\n", nangle_types);
    fprintf(otp, "\n");


    fprintf(otp, "0 %f xlo xhi\n",L[0]);
    fprintf(otp, "0 %f ylo yhi\n",L[1]);
    fprintf(otp, "0 %f zlo zhi\n",L[2]);

    fprintf(otp, "\n");

    fprintf(otp, "Masses\n\n");

	for (i = 0; i < ntypes; i++) {
        fprintf(otp, "%d %f\n");
    }
    fprintf(otp, "\n");

    fprintf(otp, "Atoms \n\n");

	for (i = 0; i < ns; i++) {
		fprintf(otp, "%d %d %d  ", i + 1, tp[i] + 1, molecID[i] + 1);
		for (j = 0; j < Dim; j++)
			fprintf(otp, "%f ", x[i][j]);

		if ( do_charges )
			fprintf(otp, "%f", charges[i]);

		for (j = Dim; j < 3; j++)
			fprintf(otp, " %f", 0.f);

		fprintf(otp, "\n");
	}
    fprintf(otp, "\n");

    if (n_total_bonds){
        fprintf(otp, "Bonds \n\n");
        int i = 0;
        for (int i1 = 0; 
        for (i = 0; i < n_total_bonds; i++) {
            fprintf(otp, "%d %d %d %d", i, 
            
        }
    }
}

*/

void read_resume( const char *nm ) {

	FILE* inp;
	inp = fopen(nm, "r");
	if (inp == NULL) {
		char death[50];
		sprintf(death, "Failed to open %s!\n", nm);
		die(death);
	}

  char tt[120];

  (void)!fgets(tt, 100, inp);
  (void)!fgets(tt, 100, inp);
  (void)!fgets(tt, 100, inp);
  int ns_tmp;
  (void)!fscanf(inp, "%d", &ns_tmp);
  if (ns != ns_tmp)
	  die("number of sites %d in resume file doens't match input!");

  (void)!fgets(tt, 100, inp);
  (void)!fgets(tt, 100, inp);
  (void)!fgets(tt, 100, inp);
  (void)!fgets(tt, 100, inp);
  (void)!fgets(tt, 100, inp);
  (void)!fgets(tt, 100, inp);


  for (int i = 0; i < ns; i++) {
	  int t1, t2, t3;
	  (void)!fscanf(inp, "%d %d %d  ", &t1, &t2, &t3);

	  for (int j = 0; j < Dim; j++)
		  (void)!fscanf(inp, "%f ", &x[i][j]);

	  for (int j = Dim; j < 3; j++) {
		  float f1;
		  (void)!fscanf(inp, "%f", &f1);
	  }
	  if (do_charges) {
		  float f1;
		  (void)!fscanf(inp, "%f", &f1);
	  }
	  
	  (void)!fgets(tt, 5, inp);
  }

  fclose(inp);
  cout << "SUCCESSFULLY READ RESUME FILE" << endl;
}
