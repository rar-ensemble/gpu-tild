#include "globals.h"
#include <sstream>

void set_defaults(void);
void read_config(const char*);
void read_charge_config(const char*);
void read_resume(const char*);
void make_group_type(int, int);
void make_group_all(void);
void write_runtime_parameters(int, string);

void read_input() {

	set_defaults();

	if (input_file.length() == 0) input_file = "input";
	ifstream in2(input_file);
	string word, line, rname;
	int config_read_flag = 0, read_resume_flag = 0; 
	int read_charge_flag = 0;
	

	while (!in2.eof()) {
		getline(in2, line);

		// Blank or commented line
		if (line.length() == 0 || line.at(0) == '#')
			continue;

		istringstream iss(line);
		// Loop over words in line
		while (iss >> word) {

			if (word == "binary_freq") {
				string toint;
				iss >> toint;
				bin_freq = stoi(toint);
			}

			else if (word == "bond") {
				if (!config_read_flag)
					die("Error in input file, bond keyword before config read");
				string toint, tofloat;
				iss >> toint;
				int btype = stoi(toint);

				string test;
				iss >> test;
				if (test != "harmonic") {
					die("You fool! That bond type is not supported!");
				}

				iss >> tofloat;
				bond_k[btype - 1] = stof(tofloat);

				iss >> tofloat;
				bond_req[btype - 1] = stof(tofloat);
			}

			else if (word == "charges") {
        
            if ( config_read_flag == 1 )
                die("charges keyword must come before read_data command!");

				read_charge_flag = 1;

				string tofloatbjerrum;
				iss >> tofloatbjerrum;
				float floatbjerrum = stof(tofloatbjerrum);
				charge_bjerrum_length = floatbjerrum;

				string tofloatlengthscale;
				iss >> tofloatlengthscale;
				float floatlengthscale = stof(tofloatlengthscale);
				charge_smearing_length = floatlengthscale;
			    do_charges = 1;
			}


			else if (word == "compute") {
				n_computes += 1;
				Computes.resize(n_computes);
				Computes[n_computes - 1].Initialize(line);
				while (!iss.eof())
					iss >> word;
			}


			else if (word == "delt") {
				string tofloat;
				iss >> tofloat;
				delt = stof(tofloat);
		        noise_mag = sqrtf(2.f * delt);
			}


			else if (word == "diffusivity") {
				if (!config_read_flag)
					die("Error in input file, diffusivities defined before config read");
				string toint, tofloat;
				iss >> toint;
				int atype = stoi(toint);
				if (atype > ntypes) {
					die("Invalid type for diffusivity!");
				}

				iss >> tofloat;
				float Df = stof(tofloat);

				Diff[atype - 1] = Df;
			}

			else if (word == "Dim") {
				string toint;
				iss >> toint;
				Dim = stoi(toint);
			}
			else if (word == "grid_freq") {
				string toint;
				iss >> toint;
				grid_freq = stoi(toint);
			}

			else if (word == "group") {
				if (config_read_flag == 0)
					die("Must read the configuration file before defining a group!");
				n_groups++;
				Groups.resize(n_groups);
				Groups.at(n_groups - 1).Initialize();
				iss >> word;
				Groups[n_groups - 1].name = word;

				iss >> word;
				Groups[n_groups - 1].command_line = line;
				if (word == "type") {
					string toint;
					iss >> toint;

					int local_type;
					local_type = stoi(toint) - 1; // shift by 1 to switch to zero indexing

					make_group_type(n_groups - 1, local_type);
				}

				else {
					char error_msg[80];
					sprintf(error_msg, "ERROR: Group type %s not supported!", word.c_str());
					die(error_msg);
				}
			}

			else if (word == "integrator") {
				//iss >> integrator;
				if (!config_read_flag)
					die("Integrator must be defined after read_data command!");
				n_integrators++;
				Integrators.resize(n_integrators);
				string grp_name, int_name;
				iss >> grp_name;
				iss >> int_name;
				Integrators[n_integrators-1].Initialize(int_name, grp_name);
			}

			else if (word == "log_freq") {
				string toint;
				iss >> toint;
				log_freq = stoi(toint);
			}
			else if (word == "max_steps") {
				string toint;
				iss >> toint;
				max_steps = stoi(toint);
			}


			else if (word == "n_erfs") {
				string toint;
				iss >> toint;
				Erfs.resize(stoi(toint));

				int i = 0;

				for (Erf& Iter : Erfs) {

					getline(in2, line);
					istringstream iss2(line);

					string convert;
					iss2 >> convert;
					if (convert != "erf") {
						cout << "Read: " << convert << " but expected 'erf' " << endl;
						die("Error in erf section of input file");
					}

					iss2 >> convert;
					Iter.type1 = stoi(convert) - 1;

					iss2 >> convert;
					Iter.type2 = stoi(convert) - 1;

					iss2 >> convert;
					Iter.initial_prefactor = stof(convert);
					Iter.final_prefactor = Iter.initial_prefactor;

					iss2 >> convert;
					Iter.Rp = stof(convert);

					iss2 >> convert;
					Iter.sigma_squared = stof(convert);
					Iter.sigma_squared *= Iter.sigma_squared;

					if (Iter.type1 >= ntypes ||
						Iter.type2 >= ntypes)
						die("Invalid particle type in Erf potential input!");

					if (iss2.tellg() != -1) {
						iss2 >> convert;
						if (convert == "ramp") {
							Iter.ramp = true;
							iss2 >> convert;
							Iter.final_prefactor = stof(convert);
							cout << "Ramping prefactor of erf style " << \
								i++ << ", between types " << Iter.type1 << " and " \
								<< Iter.type2 << " from " << Iter.initial_prefactor \
								<< " to " << Iter.final_prefactor << endl;

							cout << "Estimated per time step change: " << \
								(Iter.initial_prefactor - Iter.final_prefactor) / max_steps
								<< endl;

						}
					}
				}

			}

			else if (word == "n_fieldphases") {
				string toint;
				iss >> toint;

				Fields.resize(stoi(toint));

				int i = 0;
				for (FieldPhase& Iter : Fields) {

					getline(in2, line);
					istringstream iss2(line);

					string convert;
					iss2 >> convert;
					if (convert != "fieldphase") {
						cout << "Read: " << convert << " should be 'fieldphase'\n";
						die("Error in FieldPhase stuff in input file!\n");
					}

					iss2 >> convert;
					Iter.type1 = stoi(convert) - 1;

					iss2 >> convert;
					Iter.type2 = stoi(convert) - 1;

					iss2 >> convert;
					Iter.initial_prefactor = stof(convert);
					Iter.final_prefactor = Iter.initial_prefactor;

					iss2 >> convert;
					Iter.phase = stoi(convert);

					iss2 >> convert;
					Iter.dir = stoi(convert);

					iss2 >> convert;
					Iter.n_periods = stoi(convert);

					if (Iter.type1 >= ntypes ||
						Iter.type2 >= ntypes)
						die("Invalid particle type in FieldPhase potential input!");

					if (iss2.tellg() != -1) {
						iss2 >> convert;
						if (convert == "ramp") {
							Iter.ramp = true;
							iss2 >> convert;
							Iter.final_prefactor = stof(convert);
							cout << "Ramping prefactor of fieldphase style " << \
								i++ << ", between types " << Iter.type1 << " and " \
								<< Iter.type2 << " from " << Iter.initial_prefactor \
								<< " to " << Iter.final_prefactor << endl;

							cout << "Estimated per time step change: " << \
								(Iter.initial_prefactor - Iter.final_prefactor) / max_steps
								<< endl;

						}
					}
				}
			}

			else if (word == "n_gaussians") {
				string toint;
				iss >> toint;

				Gausses.resize(stoi(toint));


				int i = 0;
				for (Gaussian& Iter : Gausses) {

					getline(in2, line);
					istringstream iss2(line);

					string convert;
					iss2 >> convert;
					if (convert != "gaussian") {
						cout << "Read: " << convert << " should be 'gaussian'\n";
						die("Error in Gaussian stuff in input file!\n");
					}

					iss2 >> convert;
					Iter.type1 = stoi(convert) - 1;

					iss2 >> convert;
					Iter.type2 = stoi(convert) - 1;

					iss2 >> convert;
					Iter.initial_prefactor = stof(convert);
					Iter.final_prefactor = Iter.initial_prefactor;

					iss2 >> convert;
					Iter.sigma_squared = stof(convert);
					Iter.sigma_squared *= Iter.sigma_squared;

					if (Iter.type1 >= ntypes ||
						Iter.type2 >= ntypes)
						die("Invalid particle type in Gaussian input!");

					if (iss2.tellg() != -1) {
						iss2 >> convert;
						if (convert == "ramp") {
							Iter.ramp = true;
							iss2 >> convert;
							Iter.final_prefactor = stof(convert);
							cout << "Ramping prefactor of gaussian style " << \
								i++ << ", between types " << Iter.type1 << " and " \
								<< Iter.type2 << " from " << Iter.initial_prefactor \
								<< " to " << Iter.final_prefactor << endl;

							cout << "Estimated per time step change: " << \
								(Iter.initial_prefactor - Iter.final_prefactor) / max_steps
								<< endl;

						}
					}
				}
			}

			else if (word == "n_gaussian_erfs") {
				string toint;
				iss >> toint;

				GaussErfs.resize(stoi(toint));

				int i = 0;
				for (GaussianErf& Iter : GaussErfs) {

					getline(in2, line);
					istringstream iss2(line);

					string convert;
					iss2 >> convert;
					if (convert != "gaussian_erf") {
						cout << "Read: " << convert << " but expected 'gaussian_erf' " << endl;
						die("Error in gaussian_erf section of input file");
					}

					iss2 >> convert;
					Iter.type1 = stoi(convert) - 1;

					iss2 >> convert;
					Iter.type2 = stoi(convert) - 1;

					iss2 >> convert;
					Iter.initial_prefactor = stof(convert);
					Iter.final_prefactor = Iter.initial_prefactor;

					iss2 >> convert;
					Iter.Rp = stof(convert);

					iss2 >> convert;
					Iter.sigma_squared = stof(convert);
					Iter.sigma_squared *= Iter.sigma_squared;

					if (Iter.type1 >= ntypes ||
						Iter.type2 >= ntypes)
						die("Invalid particle type in GaussianErf potential input!");

					if (iss2.tellg() != -1) {
						iss2 >> convert;
						if (convert == "ramp") {
							Iter.ramp = true;
							iss2 >> convert;
							Iter.final_prefactor = stof(convert);
							cout << "Ramping prefactor of gaussian_erf style " << \
								i++ << ", between types " << Iter.type1 << " and " \
								<< Iter.type2 << " from " << Iter.initial_prefactor \
								<< " to " << Iter.final_prefactor << endl;

							cout << "Estimated per time step change: " << \
								(Iter.initial_prefactor - Iter.final_prefactor) / max_steps
								<< endl;

						}
					}
				}

			}




			else if (word == "Nx") {
				string toint;
				iss >> toint;
				Nx[0] = stoi(toint);
			}

			else if (word == "Ny") {
				string toint;
				iss >> toint;
				Nx[1] = stoi(toint);
			}

			else if (word == "Nz") {
				string toint;
				iss >> toint;
				if (Dim > 2)
					Nx[2] = stoi(toint);
			}



			else if (word == "pmeorder") {
				string toint;
				iss >> toint;
				pmeorder = stoi(toint);
				if (pmeorder > 4)
					die("PMEORDER greater than 4 not supported!");
			}

			else if (word == "rand_seed") {
				string toint;
				iss >> toint;
				RAND_SEED = stoi(toint);
			}


			else if (word == "read_data") {
				string to_char;
				iss >> to_char;
				if (read_charge_flag == 1) {
					cout << "READING CHARGE CONFIG" << endl;
					read_charge_config(to_char.c_str());
				}
				else {
					read_config(to_char.c_str());
				}
				config_read_flag = 1;

				Groups.resize(1);
				Groups.at(0).Initialize();

				// Initialize the group "all"
				make_group_all();
			}
			else if (word == "read_resume") {
				string to_char;
				iss >> to_char;
				read_resume(to_char.c_str());
				read_resume_flag = 1;
				rname = to_char;
			}




			else if (word == "skip_steps") {
				string toint;
				iss >> toint;
				skip_steps = stoi(toint);
			}

			else if (word == "struc_freq") {
				string toint;
				iss >> toint;
				struc_freq = stoi(toint);
			}

			else if (word == "threads") {
				string toint;
				iss >> toint;
				threads = stoi(toint);
			}

			else if (word == "traj_freq") {
				string toint;
				iss >> toint;
				traj_freq = stoi(toint);
			}

			else if (word == "traj_name") {
				dump_name.clear();
				iss >> dump_name;
			}



			else {
				cout << "Invalid keyword: " << word << endl;
				die("Invalid keyword error");
			}
		}// while (iss >> word)
		cout << line << endl;
	}// while (!in2.eof())

  if ( using_GJF && delt < 0.002f ) {
    cout << "WARNING!!!!! WARNING!!!! WARNING!!!!" << endl;
    cout << "Previous tests have seen problems with GJF algorithm and time steps < 0.002" << endl;
  }


	write_runtime_parameters(read_resume_flag, rname);
	
}


void set_defaults() {
	
	Dim = 2;
	delt = 0.005f;
	pmeorder = 1;
	for (int j = 0; j < Dim; j++) {
		Nx[j] = 35;
	}

	threads = 512;
	noise_mag = sqrtf(2.0f * delt);
	n_groups = 1;
  n_integrators = 0;
	n_computes = 0;
	using_GJF = 0;
  do_charges = 0;

	max_steps = 100000;
	RAND_SEED = int(time(0));
	log_freq = 2500;
	grid_freq = 0;
	traj_freq = 0;
	struc_freq = 0;
	skip_steps = 20000;
	bin_freq = 5000;
	dump_name = "traj.lammpstrj";
}


void write_runtime_parameters(int rflag, string rname) {
		ofstream otp;
		otp.open("NEED_A_NAME.input");

		otp << "Dim " << Dim << endl;
		otp << "Nx " << Nx[0] << endl;
		otp << "Ny " << Nx[1] << endl;
		if (Dim == 3) {
			otp << "Nz " << Nx[2] << endl;
		}
		
		otp << "delt " << delt << endl;

		for (int i = 1; i < n_groups; i++) {
			otp << Groups[i].command_line << endl;
		}

		for (int i = 0; i < n_integrators; i++) {
			int grp_id = Integrators.at(i).group_index;
			otp << "integrator " << Groups.at(grp_id).name << " " << Integrators.at(i).name << endl;
		}

        otp << "rand_seed " << RAND_SEED << endl;
		if (rflag)
			otp << "read_resume " << rname << endl;
		otp << "max_steps " << max_steps << endl;
		otp << "log_freq " << log_freq << endl;
		otp << "threads " << threads << endl;
		otp << "grid_freq " << grid_freq << endl;
		otp << "traj_freq " << traj_freq << endl;
		otp << "binary_freq " << bin_freq << endl;
		otp << "struc_freq " << struc_freq << endl;
		otp << "skip_steps " << skip_steps << endl;
		otp << "pmeorder " << pmeorder << endl;
		if (do_charges == 1)
			otp << "charges " << charge_bjerrum_length << " " << charge_smearing_length << endl;

		for (int i = 0; i < ntypes; i++) {
			otp << "diffusivity " << i + 1 << " " << Diff[i] << endl;
		}
		
		
		if (!Gausses.empty()) {
		otp << "n_gaussians " << Gausses.size() << endl;
			for (Gaussian& Iter : Gausses) {
				otp << "gaussian" << Iter.type1 + 1 << " " \
					<< Iter.type2 + 1 << " " << Iter.initial_prefactor \
					<< " " << Iter.sigma_squared<< endl;
			}

		}

		if (!Fields.empty()) {
			otp << "n_fieldphases " << Fields.size() << endl;
			for (FieldPhase& Iter : Fields) {
				otp << "fieldphase " << Iter.type1 + 1 << " " \
					<< Iter.type2 + 1 << " " << Iter.initial_prefactor \
					<< " " << Iter.phase << " " << Iter.dir  \
					<< " " << Iter.n_periods << endl;
			}
		}

		if (!Erfs.empty()) {
			otp << "n_erfs " << Erfs.size() << endl;
			for (Erf& Iter : Erfs) {
				otp << "erf " << Iter.type1 + 1 << " " << Iter.type2 + 1 << " " \
					<< Iter.initial_prefactor << " " << Iter.Rp << " " << Iter.sigma_squared << endl;
			}
		otp << endl;
		}

		if (!GaussErfs.empty()) {
			otp << "n_gaussian_erfs " << GaussErfs.size() << endl;
			for (GaussianErf& Iter : GaussErfs) {
				otp << "erf " << Iter.type1 + 1 << " " << Iter.type2 + 1 << " " \
					<< Iter.initial_prefactor << " " << Iter.Rp << " " << Iter.sigma_squared << endl;
			}
			otp << endl;
		}

		for (int i = 0; i < n_computes; i++)
			otp << Computes.at(i).input_command << endl;

		for (int i = 0; i < nbond_types; i++)
			otp << "bond " << i + 1 << " harmonic " << bond_k[i] \
			<< " " << bond_req[i] << endl;
}
