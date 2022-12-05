// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2021           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo 
// Corey Oses
// Frisco Rose (AFLUX)
// Marco Esters (AFLOW DB)

#ifndef _AFLOWLIB_H_
#define _AFLOWLIB_H_

#include "aflow.h"
#include "SQLITE/aflow_sqlite.h"
//[OBSOLETE] [KY] #include "aflow_contrib_kesong.h" //CO20180515

using std::vector;
using std::string;

#define NOSG string("NNN #0")
#define CHMODWEB FALSE
#define INF 1E9
#define ENERGY_ATOM_ERROR_meV 50
#define PRESSURE_ZERO_ENTHALPY_ENERGY _FLOAT_TOL_ //1e-6

extern vector<string> vLibrary_LIB2;   // need to remove somwhoe
extern vector<string> vLibrary_LIB3;
extern vector<string> vLibrary_LIB1;
extern vector<string> vLibrary_ICSD;
extern vector<string> vLibrary_AURO;
extern vector<vector<string> > vLibrary_LIB2_tokens;
extern vector<vector<string> > vLibrary_LIB3_tokens;
extern vector<vector<string> > vLibrary_LIB1_tokens;
extern vector<vector<string> > vLibrary_ICSD_tokens;
extern vector<vector<string> > vLibrary_AURO_tokens;

// namespace for the web
// ***************************************************************************
namespace aflowlib {
  class _aflowlib_entry {
    public:
      // constructor destructor                                 // constructor/destructor
      _aflowlib_entry();                                        // default, just allocate
      ~_aflowlib_entry();                                       // kill everything
      _aflowlib_entry(const _aflowlib_entry& b);                // constructor copy
      _aflowlib_entry(const string& file);                      // constructor from file
      const _aflowlib_entry& operator=(const _aflowlib_entry &b); // copy
      // CONTROL
      string entry;vector<string> ventry;                       // ventry split by "|"
      string auid;deque<string> vauid;                          // AFLOW UNIQUE IDENTIFIER // AFLOW UNIQUE IDENTIFIER SPLIT
      string aurl;deque<string> vaurl;                          // AFLOW RESEARCH LOCATOR and TOKENS
      string system_name;                                       //ME20190125 - system_name of the calculation
      string keywords;deque<string> vkeywords;                  // keywords inside
      string aflowlib_date;vector<string> vaflowlib_date;       // CONTAINS 2 DATES: [0]=lock_date, [1]=lib2raw_date
      string aflowlib_version;                                  // version
      string aflowlib_entries;vector<string> vaflowlib_entries; // this contains the subdirectories that can be associated
      int aflowlib_entries_number;                              // their number
      string aflow_version;                                     // version
      string catalog;                                           // ICSD,LIB2, etc.
      string data_api,data_source,data_language;                // version/source/language
      string error_status;                                            // ERROR ??
      string author;vector<string> vauthor; 
      int calculation_cores;double calculation_memory,calculation_time;
      string corresponding;vector<string> vcorresponding; 
      string loop;vector<string> vloop;                         // postprocessing
      int node_CPU_Cores;                                       // computer
      double node_CPU_MHz;                                      // computer
      string node_CPU_Model;                                    // computer
      double node_RAM_GB;                                       // computer
      // materials
      string Bravais_lattice_orig,Bravais_lattice_relax;        // structure
      string code;
      string composition;vector<double> vcomposition;
      string compound;
      double density;
      double density_orig; //DX20190124 - add original crystal info
      string dft_type;vector<string> vdft_type;
      double eentropy_cell,eentropy_atom;
      double Egap,Egap_fit;
      string Egap_type;  
      double energy_cell,energy_atom,energy_atom_relax1;
      double energy_cutoff;
      double delta_electronic_energy_convergence;
      double delta_electronic_energy_threshold;
      uint nkpoints,nkpoints_irreducible,kppra;
      string kpoints;
      xvector<int> kpoints_nnn_relax,kpoints_nnn_static;
      vector<string> kpoints_pairs;
      double kpoints_bands_path_grid;
      double enthalpy_cell,enthalpy_atom;
      double enthalpy_formation_cell,enthalpy_formation_atom;
      double enthalpy_formation_cce_300K_cell,enthalpy_formation_cce_300K_atom; //CO20200624
      double enthalpy_formation_cce_0K_cell,enthalpy_formation_cce_0K_atom; //CO20200624
      double entropy_forming_ability; //CO20200624
      double entropic_temperature;
      string files;vector<string> vfiles;
      string files_LIB;vector<string> vfiles_LIB;
      string files_RAW;vector<string> vfiles_RAW;
      string files_WEB;vector<string> vfiles_WEB;
      string forces;vector<xvector<double> > vforces;
      string geometry;vector<double> vgeometry; // a,b,c and unit_cell_angles (b,c) (a,c) (a,b)
      string geometry_orig;vector<double> vgeometry_orig; // a,b,c and unit_cell_angles (b,c) (a,c) (a,b) //DX20190124 - add original crystal info
      string lattice_system_orig,lattice_variation_orig,lattice_system_relax,lattice_variation_relax;
      string ldau_TLUJ;
      vector<vector<double> > vLDAU;  //ME20190124
      uint natoms;
      uint natoms_orig; //DX20190124 - add original crystal info
      string nbondxx;vector<double> vnbondxx;
      uint nspecies;
      string Pearson_symbol_orig,Pearson_symbol_relax;
      string positions_cartesian;vector<xvector<double> > vpositions_cartesian;
      string positions_fractional;vector<xvector<double> > vpositions_fractional;
      double pressure;           // the true applied pressure (PSTRESS)
      string stress_tensor;vector<double> vstress_tensor; // (1,1),(1,2),(1,3),(2,1),(2,2),(2,3),(3,1),(3,2),(3,3)
      double pressure_residual;  // the leftover pressure due to convergence 
      double Pulay_stress;       // the leftover pressure due to incomplete basis set
      string prototype;
      double PV_cell,PV_atom;
      double scintillation_attenuation_length;
      string sg,sg2;vector<string> vsg,vsg2; //CO20180101
      uint spacegroup_orig,spacegroup_relax;
      string species;vector<string> vspecies;
      string species_pp;vector<string> vspecies_pp;
      string species_pp_version;vector<string> vspecies_pp_version;
      string species_pp_ZVAL;vector<double> vspecies_pp_ZVAL;
      string species_pp_AUID;vector<string> vspecies_pp_AUID;
      string METAGGA; // empty if none, potential type is in xPOTCAR/xOUTCAR
      double spin_cell,spin_atom;
      string spinD;vector<double> vspinD;
      string spinD_magmom_orig;vector<double> vspinD_magmom_orig;
      double spinF;
      string sponsor;vector<string> vsponsor; 
      string stoichiometry;vector<double> vstoichiometry;
      double valence_cell_std,valence_cell_iupac;
      double volume_cell,volume_atom;
      double volume_cell_orig,volume_atom_orig; //DX20190124 - add original crystal info
      //DX20190124 - added original symmetry info - START
      // SYMMETRY
      string crystal_family_orig;
      string crystal_system_orig;
      string crystal_class_orig;
      string point_group_Hermann_Mauguin_orig;
      string point_group_Schoenflies_orig;
      string point_group_orbifold_orig;
      string point_group_type_orig;
      uint point_group_order_orig;
      string point_group_structure_orig;
      string Bravais_lattice_lattice_type_orig;
      string Bravais_lattice_lattice_variation_type_orig;
      string Bravais_lattice_lattice_system_orig;
      string Bravais_superlattice_lattice_type_orig;
      string Bravais_superlattice_lattice_variation_type_orig;
      string Bravais_superlattice_lattice_system_orig;
      string Pearson_symbol_superlattice_orig;
      string reciprocal_geometry_orig;vector<double> vreciprocal_geometry_orig;
      double reciprocal_volume_cell_orig;
      string reciprocal_lattice_type_orig;
      string reciprocal_lattice_variation_type_orig;
      string Wyckoff_letters_orig;
      string Wyckoff_multiplicities_orig;
      string Wyckoff_site_symmetries_orig;
      //DX20190124 - added original symmetry info - END
      //DX20180823 - added more symmetry info - START
      // SYMMETRY
      string crystal_family;
      string crystal_system;
      string crystal_class;
      string point_group_Hermann_Mauguin;
      string point_group_Schoenflies;
      string point_group_orbifold;
      string point_group_type;
      uint point_group_order;
      string point_group_structure;
      string Bravais_lattice_lattice_type;
      string Bravais_lattice_lattice_variation_type;
      string Bravais_lattice_lattice_system;
      string Bravais_superlattice_lattice_type;
      string Bravais_superlattice_lattice_variation_type;
      string Bravais_superlattice_lattice_system;
      string Pearson_symbol_superlattice;
      string reciprocal_geometry;vector<double> vreciprocal_geometry;
      double reciprocal_volume_cell;
      string reciprocal_lattice_type;
      string reciprocal_lattice_variation_type;
      string Wyckoff_letters;
      string Wyckoff_multiplicities;
      string Wyckoff_site_symmetries;
      //DX20180823 - added more symmetry info - END
      //DX20190208 - added anrl info - START
      string aflow_prototype_label_orig; //DX20201001 - renamed
      string aflow_prototype_params_list_orig; //DX20201001 - renamed
      string aflow_prototype_params_values_orig; //DX20201001 - renamed
      string aflow_prototype_label_relax; //DX20201001 - renamed
      string aflow_prototype_params_list_relax; //DX20201001 - renamed
      string aflow_prototype_params_values_relax; //DX20201001 - renamed
      //DX20190208 - added anrl info - END
      string pocc_parameters; //CO20200731
      // AGL/AEL
      double agl_thermal_conductivity_300K;//  (W/m*K)
      double agl_debye;//  (K)
      double agl_acoustic_debye;//  (K)
      double agl_gruneisen;// 
      double agl_heat_capacity_Cv_300K;//  (kB/cell)
      double agl_heat_capacity_Cp_300K;//  (kB/cell)
      double agl_thermal_expansion_300K;//  (1/K)
      double agl_bulk_modulus_static_300K;//  (GPa)
      double agl_bulk_modulus_isothermal_300K;//  (GPa)
      string agl_poisson_ratio_source;//
      double agl_vibrational_free_energy_300K_cell;// (meV/cell) //CT20181212
      double agl_vibrational_free_energy_300K_atom;// (meV/atom) //CT20181212
      double agl_vibrational_entropy_300K_cell;// (meV/cell*K) //CT20181212
      double agl_vibrational_entropy_300K_atom;// (meV/atom*K) //CT20181212
      double ael_poisson_ratio ;//
      double ael_bulk_modulus_voigt;//  (GPa)
      double ael_bulk_modulus_reuss;//  (GPa)
      double ael_shear_modulus_voigt;//  (GPa)
      double ael_shear_modulus_reuss;//  (GPa)
      double ael_bulk_modulus_vrh;//  (GPa)
      double ael_shear_modulus_vrh;//  (GPa)
      double ael_elastic_anisotropy;// //CO20181128
      double ael_youngs_modulus_vrh;//  (GPa) //CT20181212
      double ael_speed_sound_transverse;// (m/s) //CT20181212
      double ael_speed_sound_longitudinal;// (m/s) //CT20181212
      double ael_speed_sound_average;// (m/s) //CT20181212
      double ael_pughs_modulus_ratio;// //CT20181212
      double ael_debye_temperature;// (K) //CT20181212
      double ael_applied_pressure;// (GPa) //CT20181212
      double ael_average_external_pressure; // (GPa) //CT20181212
      xmatrix<double> ael_stiffness_tensor;  //ME20191105
      xmatrix<double> ael_compliance_tensor;  //ME20191105
      // QHA  //AS20200831
      double gruneisen_qha; //AS20200831
      double gruneisen_qha_300K; //AS20200903
      double thermal_expansion_qha_300K; //AS20200831
      double modulus_bulk_qha_300K; //AS20200831
      double modulus_bulk_derivative_pressure_qha_300K; //AS20201008
      double heat_capacity_Cv_atom_qha_300K; //AS20201008
      double heat_capacity_Cv_cell_qha_300K; //AS20201207
      double heat_capacity_Cp_atom_qha_300K; //AS20201008
      double heat_capacity_Cp_cell_qha_300K; //AS20201207
      double volume_atom_qha_300K; //AS20201008
      double energy_free_atom_qha_300K; //AS20201008
      double energy_free_cell_qha_300K; //AS20201207
      // BADER
      string bader_net_charges;vector<double> vbader_net_charges;//electrons
      string bader_atomic_volumes;vector<double> vbader_atomic_volumes;//Angst^3
      // legacy
      string server;vector<string> vserver;vector<vector<string> > vserverdir; 
      string icsd;
      string stoich;vector<double> vstoich;
      // apennsy
      string structure_name;                            // name for apennsy
      string structure_description;                     // description for apennsy
      double distance_gnd;                              // distance_gnd
      double distance_tie;                              // distance_tie
      bool pureA,pureB;                                 // pureA,pureB
      bool fcc,bcc,hcp;                                 // options for lattices
      double stoich_a,stoich_b;                         // stoich_b,stoich_b
      double bond_aa,bond_ab,bond_bb;                   // bond_xx // BOND_XX [norm V_ATOM^0.33]
      vector<uint> vNsgroup;                            // vNsgroups
      vector<string> vsgroup;                           // vsgroups
      vector<xstructure> vstr;                          // vstructures
      // functions
      bool FixDescription(void);                        // fix description names
      void GetSGROUP(string aflowlibentry);             // disassemble SG
      void clear();                                              // free space
      uint Load(const stringstream& stream,ostream& oss);        // load from stringstream it std is cout
      uint Load(const string& entry,ostream& oss);               // load from string it std is cout
      uint file2aflowlib(const string& file,ostream& oss=std::cout);       // load from file
      uint url2aflowlib(const string& url,ostream& oss,bool=TRUE); // load from the web (VERBOSE)
      string aflowlib2string(string="out", bool=false);          //ME20210408 - added PRINT_NULL
      string aflowlib2file(string file,string="out");            //
      string POCCdirectory2MetadataAUIDjsonfile(const string& directory,uint salt=0);           //CO20200624 - get contents of auid_metadata.json 
      bool directory2auid(const string& directory);                                         // from directory and AURL gives AUID and VAUID
      double enthalpyFormationCell(int T=300) const;                                        //CO20200624 - CCE correction added
      double enthalpyFormationAtom(int T=300) const;                                        //CO20200624 - CCE correction added
      void correctBadDatabase(bool verbose=true,ostream& oss=cout);                         //CO20171202 - apennsy fixes
      void correctBadDatabase(ofstream& FileMESSAGE,bool verbose=true,ostream& oss=cout);   //CO20171202 - apennsy fixes
      bool ignoreBadDatabase() const;                                                       //CO20171202 - apennsy fixes
      bool ignoreBadDatabase(string& reason) const;                                         //CO20171202 - apennsy fixes
      string getPathAURL(ostream& oss=cout, bool load_from_common=false);                   // converts entry.aurl to url/path (common)
      string getPathAURL(ofstream& FileMESSAGE, ostream& oss, bool load_from_common=false); // converts entry.aurl to url/path (common)
      vector<string> getSpeciesAURL(ostream& oss);                                          //CO20210201 - extracts species from aurl
      vector<string> getSpeciesAURL(ofstream& FileMESSAGE,ostream& oss);                    //CO20210201 - extracts species from aurl
      //ML stoich features
      void getStoichFeatures(vector<string>& vheaders,const string& e_props=_AFLOW_XELEMENT_PROPERTIES_ALL_);
      void getStoichFeatures(vector<string>& vheaders,vector<double>& vfeatures,bool vheaders_only=false,const string& e_props=_AFLOW_XELEMENT_PROPERTIES_ALL_);
    private:                                                     //
      void free();                                               // free space
      void copy(const _aflowlib_entry& b);                       //
  };
}

namespace aflowlib {
  string VASPdirectory2auid(const string& directory,const string& aurl);   //CO20200624 - moving from inside _aflowlib_entry
  uint auid2vauid(const string auid, deque<string>& vauid);                // splits the auid into vauid
  string auid2directory(const string auid);                                // gives AUID directory from existence of vauid
  void insertStoichStats(const vector<string> vstats,const xvector<double>& nspecies_xv,const xvector<double>& stoich_xv,vector<double>& vfeatures);  //CO20201111
}

// ***************************************************************************
// AFLUX STUFF
namespace aflowlib {
  class APIget {
    private:
      struct sockaddr_in client;
      int sock;
      int PORT;// = 80; //CO20180401
      string Summons;
      string API_Path;
      string Domain;
      bool establish();
    public:
      //DX20210615 [OBSOLETE] APIget( string a_Summons="", string a_API_Path="/search/API/?", string a_Domain="aflowlib.duke.edu" ): PORT(80), Summons(a_Summons), API_Path(a_API_Path), Domain(a_Domain) {}; //CO20181226
      APIget( string a_Summons="", string a_API_Path="/API/aflux/?", string a_Domain="aflow.org" ): PORT(80), Summons(a_Summons), API_Path(a_API_Path), Domain(a_Domain) {}; //CO20181226 //DX20210615 - updated domain name
      void reset( string a_Summons="#", string a_API_Path="", string a_Domain="" );
      friend ostream& operator<<( ostream& output, APIget& a );
  };
}

// ***************************************************************************
namespace aflowlib {
  bool json2aflowlib(const string& json,string key,string& value); //ME I do not understand why you have out this function as a private member... it would have been simpler to be able to use everywhere
  uint auid2present(string auid,string& aurl,int mode=1); // returns json.size() if found...
  bool AflowlibLocator(const string& in,string& out,const string& mode);
  string AflowlibLocator(string options,string mode);
  string AFLUXCall(const aurostd::xoption& vpflow); //DX20190206 - add AFLUX functionality for command line 
  string AFLUXCall(const vector<string>& matchbook); //DX20190206 - add AFLUX functionality
  string AFLUXCall(const string& summons); //DX20190206 - add AFLUX functionality 
  vector<vector<std::pair<string,string> > > getPropertiesFromAFLUXResponse(const string& response); //DX20190206 - get properties from AFLUX response
  string getSpaceGroupAFLUXSummons(const vector<uint>& space_groups, uint relaxation_step); //DX20200929
  string getSpaceGroupAFLUXSummons(uint space_group_number, uint relaxation_step, bool only_one_sg=true); //DX20200929
  // [OBSOLETE] uint WEB_Aflowlib_Entry_PHP(string options,ostream& oss); //SC20200327
  uint WEB_Aflowlib_Entry(string options,ostream& oss); 
  // [OBSOLETE] uint WEB_Aflowlib_Entry_PHP3(string options,ostream& oss);  //SC20190813
}

// ***************************************************************************
// to create/analyze the tokens and load the stuff up for qhull
namespace aflowlib {
  bool TokenPresentAFLOWLIB(string line,string query);
  string TokenExtractAFLOWLIB(string line,string query);
  bool TokenExtractAFLOWLIB(string line,string query,string &value);
  bool TokenExtractAFLOWLIB(string line,string query,int &value);
  bool TokenExtractAFLOWLIB(string line,string query,uint &value);
  bool TokenExtractAFLOWLIB(string line,string query,float &value);
  bool TokenExtractAFLOWLIB(string line,string query,double &value);
  bool TokenExtractAFLOWLIB(string line,string query,vector<string> &value);
  bool TokenExtractAFLOWLIB(string line,string query,vector<int> &value);
  bool TokenExtractAFLOWLIB(string line,string query,vector<uint> &value);
  bool TokenExtractAFLOWLIB(string line,string query,vector<float> &value);
  bool TokenExtractAFLOWLIB(string line,string query,vector<double> &value);
}

namespace aflowlib {
  uint LOAD_Library_LIBRARY(string file,string="",bool=FALSE);
  uint LOAD_Library_ALL(string options,string="",bool=FALSE);
  uint LOAD_Library_ALL(string options,bool=FALSE);
  uint LOAD_Library_ALL(bool=FALSE);
}

namespace aflowlib {
  uint GREP_Species_ALL(vector<string> vspecies,                         // IN  [0,nspecies[ the species Ag,Cd,...   nspecies=number of these items    nspecies=naries
      vector<string>& vspecies_pp,                     // IN  [0,nspecies[ the pseudopotentials Ag_pv, Cd_sv
      vector<vector<string> >& vList,                  // OUT [0,naries[*[0,vList.size()[ returns the lines of the library containing A,B,C,AB,AC,BC,ABC....
      vector<vector<vector<string> > > &vList_tokens,  // OUT [0,naries[*[0,vList.size()[*[0,size_tokens[ returns the tokens for each line of the vList
      vector<vector<vector<string> > > &vList_species, // OUT [0,naries[*[0,vList.size()[*nspecies  returns the species present
      vector<double> &vList_Hmin,                      // OUT [0,nspecies[ returns the min enthalpy for reference
      vector<string> &vList_Pmin,                      // OUT [0,nspecies[ returns the prototype for reference
      vector<uint> &vList_Imin,                        // OUT [0,nspecies[ returns the line index of vList.at(ispecies), in which we have the min enthalpy for reference
      vector<vector<vector<double> > > &vList_concs,   // OUT [0,naries[*[0,vList.size()[*[0.nspecies[ the concentrations AxAyCz... where x+y+z=1 and it contains also ZEROS so that 0 0.25 0.75 is allowed
      vector<vector<double> > &vList_Ef);              // OUT [0,naries[*[0,vList.size()[ returns the formation energy of the list
  uint GREP_Species_ALL(string species,                                  // IN  the species Ag,...   nspecies=number of these items
      string& species_pp,                              // IN  the pseudopotentials Ag_pv,
      vector<string>& vList,                           // OUT [0,vList.size()[ returns the lines of the library containing A,B,C,AB,AC,BC,ABC....
      vector<vector<string> > &vList_tokens,           // OUT [0,vList.size()[*[0,size_tokens[ returns the tokens for each line of the vList
      double& List_Hmin,                               // OUT returns the min enthalpy for reference
      string& List_Pmin,                               // OUT returns the prototype for reference
      uint& List_Imin);                                // OUT returns the line index of vList, in which we have the min enthalpy for reference
  uint GREP_Species_ALL(vector<string> vspecies,vector<string>& vspecies_pp,vector<vector<string> >& vList,vector<vector<vector<string> > > &vList_tokens);
  uint GREP_Species_ALL(vector<string> vspecies,vector<string>& vspecies_pp,vector<vector<string> >& vList);
  uint GREP_Species_ALL(string species,                                  // IN  the species Ag,...   nspecies=number of these items
      string& species_pp,                              // IN  the pseudopotentials Ag_pv,
      double& List_Hmin,                               // OUT returns the min enthalpy for reference
      string& List_Pmin);                              // OUT returns the prototype for reference
}

namespace aflowlib {
  bool ALIBRARIES(string options);
  bool LIBS_EFormation(string options);
}

// ***************************************************************************
// LIBS to RAWS of each entry
namespace aflowlib {
  uint GetSpeciesDirectory(const string& directory,vector<string>& vspecies);
}
namespace aflowlib {
  void XFIX_LIBRARY_ALL(string LIBRARY_IN,vector<string>);
  // void LIB2RAW_LIBRARY_ALL(string LIBRARY_IN);
  string LIB2RAW_CheckProjectFromDirectory(const string& directory);
  bool LIB2RAW_ALL(string options,bool overwrite);
  bool LIB2RAW_FileNeeded(string directory_LIB,string fileLIB,string directory_RAW,string fileRAW,vector<string> &vfiles,const string& MESSAGE);
  // [OBSOLETE] bool LIB2RAW(vector<string> argv,bool overwrite);
  void CleanDirectoryLIB(string& directory);  //CO20200624
  bool LIB2RAW(const string& options,bool overwrite,bool LOCAL=false);
  bool XPLUG_CHECK_ONLY(const vector<string>& argv); //CO20200501
  bool XPLUG_CHECK_ONLY(const vector<string>& argv,deque<string>& vdirsOUT,deque<string>& vzips,deque<string>& vcleans); //CO20200501
  bool XPLUG(const vector<string>& argv);  //CO20200501
  bool AddFileNameBeforeExtension(string _file,string addendum,string& out_file); //CO20171025
  bool LIB2RAW_Calculate_FormationEnthalpy(aflowlib::_aflowlib_entry& data,const xstructure& xstr,const string& MESSAGE); //CO20200731
  bool LIB2RAW_Loop_Thermodynamics(const string& directory_LIB,const string& directory_RAW,vector<string> &vfiles,aflowlib::_aflowlib_entry&,const string& MESSAGE,bool LOCAL=false);
  // [OBSOLETE]  bool LIB2RAW_Loop_DATA(const string& directory_LIB,const string& directory_RAW,vector<string> &vfiles,aflowlib::_aflowlib_entry& data,const string& MESSAGE);
  bool LIB2RAW_Loop_Static(const string& directory_LIB,const string& directory_RAW,vector<string> &vfiles,aflowlib::_aflowlib_entry&,const string& MESSAGE);  //CO20200731
  bool LIB2RAW_Loop_Bands(const string& directory_LIB,const string& directory_RAW,vector<string> &vfiles,aflowlib::_aflowlib_entry&,const string& MESSAGE);
  bool LIB2RAW_Loop_Magnetic(const string& directory_LIB,const string& directory_RAW,vector<string> &vfiles,aflowlib::_aflowlib_entry&,const string& MESSAGE);
  bool LIB2RAW_Loop_Bader(const string& directory_LIB,const string& directory_RAW,vector<string> &vfiles,aflowlib::_aflowlib_entry&,const string& MESSAGE);
  bool LIB2RAW_Loop_AGL(const string& directory_LIB,const string& directory_RAW,vector<string> &vfiles,aflowlib::_aflowlib_entry&,const string& MESSAGE);
  bool LIB2RAW_Loop_AEL(const string& directory_LIB,const string& directory_RAW,vector<string> &vfiles,aflowlib::_aflowlib_entry&,const string& MESSAGE);
  bool LIB2RAW_Loop_QHA(const string& directory_LIB,const string& directory_RAW,vector<string> &vfiles,aflowlib::_aflowlib_entry& data,const string& MESSAGE);  //AS20200831
  bool LIB2RAW_Loop_LOCK(const string& directory_LIB,const string& directory_RAW,vector<string> &vfiles,aflowlib::_aflowlib_entry& data,const string& MESSAGE);
  bool LIB2RAW_Loop_POCC(const string& directory_LIB,const string& directory_RAW,vector<string> &vfiles,aflowlib::_aflowlib_entry& data,const string& MESSAGE);  //CO20200624
  bool LIB2RAW_Loop_PATCH(const string& directory_LIB,const string& directory_RAW,vector<string> &vfiles,aflowlib::_aflowlib_entry& data,const string& MESSAGE);
  bool LIB2LIB(const string& options,bool overwrite,bool LOCAL=false); //CT20181212
}

namespace aflowlib {
  bool VaspFileExist(const string& str_dir, const string& FILE);
  string vaspfile2stringstream(const string& str_dir, const string& FILE, stringstream& streamFILE);
  string vaspfile2stringstream(const string& str_dir, const string& FILE);
}

// ***************************************************************************
#define HTRESOURCE_MODE_NONE  0
#define HTRESOURCE_MODE_PHP_AUTHOR   4
#define HTRESOURCE_MODE_PHP_THRUST   5
#define HTRESOURCE_MODE_PHP_ALLOY    6

// ***************************************************************************
// _OUTREACH CLASS
class _outreach {
  public:
    // constructors/destructors
    _outreach(void);    // do nothing
    _outreach(const _outreach& b);    // do nothing
    ~_outreach();        // do nothing
    // OPERATORS                                              // --------------------------------------
    const _outreach& operator=(const _outreach& b);             // some operators
    void clear(void);                                         // clear
    // CONTENT
    // [OBSOLETE] string print_mode; // "TXT","LATEX","HTML"
    // [OBSOLETE] inside  XHOST.vflag_control.flag("PRINT_MODE::TXT");
    // [OBSOLETE] inside  XHOST.vflag_control.flag("PRINT_MODE::LATEX");
    // [OBSOLETE] inside  XHOST.vflag_control.flag("PRINT_MODE::HTML");
    // [OBSOLETE] bool print_doi; inside  XHOST.vflag_control.flag("PRINT_MODE::DOI");
    // [OBSOLETE] bool print_bibtex; inside  XHOST.vflag_control.flag("PRINT_MODE::BIBTEX"); //SC20201228
    // [OBSOLETE] bool print_pdf; inside  XHOST.vflag_control.flag("PRINT_MODE::PDF");
    // [OBSOLETE] bool print_wnumber; inside  XHOST.vflag_control.flag("PRINT_MODE::NUMBER");
    uint wnumber;
    bool newflag;
    uint year;
    vector<string> vauthor;
    string title;
    string journal,link,arxiv,supplementary; 
    string place,date;
    string type;   // ARTICLE PRESENTATION_TALK PRESENTATION_SEMINAR PRESENTATION_COLLOQUIUM PRESENTATION_KEYNOTE PRESENTATION_PLENARY PRESENTATION_TUTORIAL PRESENTATION_CONTRIBUTED PRESENTATION_POSTER
    bool _isinvited;       // YES
    bool _isonline;       // YES
    string host;           // who invited
    string abstract;       // if available and in LaTeX
    string pdf;
    string doi;
    string bibtex,bibtex_journal,bibtex_volume,bibtex_issue,bibtex_pages,bibtex_year;  //SC20201228
    vector<string> vextra_html,vextra_latex;
    vector<string> vkeyword,vsponsor,valloy;
    // operators/functions                                    // operator/functions
    friend ostream& operator<<(ostream &,const _outreach&);   // print
    string print_string(void);                                // print string from <<
  private:                                                    // ---------------------------------------
    void free();                                              // to free everything
    void copy(const _outreach& b);                            // the flag is necessary because sometimes you need to allocate the space.
};

uint voutreach_load(vector<_outreach>& voutreach,string what2print);
void voutreach_print(uint _mode,ostream& oss,string what2print);
void voutreach_print_everything(ostream& oss,const vector<string>& vitems,string msg1,string msg2,string sectionlabel);

// sort

class _sort_outreach_outreach_year_ {
  public:
    bool operator()(const _outreach& a1, const _outreach& a2) const {
      return (bool) (a1.year<a2.year);}}; // sorting through reference
class _rsort_outreach_outreach_year_ {
  public:
    bool operator()(const _outreach& a1, const _outreach& a2) const {
      return (bool) (a1.year>a2.year);}}; // sorting through reference
class _sort_outreach_outreach_wnumber_ {
  public:
    bool operator()(const _outreach& a1, const _outreach& a2) const {
      return (bool) (a1.wnumber<a2.wnumber);}}; // sorting through reference
class _rsort_outreach_outreach_wnumber_ {
  public:
    bool operator()(const _outreach& a1, const _outreach& a2) const {
      return (bool) (a1.wnumber>a2.wnumber);}}; // sorting through reference

uint voutreach_sort_year(vector<_outreach>& voutreach);
uint voutreach_rsort_year(vector<_outreach>& voutreach);
uint voutreach_sort_wnumber(vector<_outreach>& voutreach);
uint voutreach_rsort_wnumber(vector<_outreach>& voutreach);

uint voutreach_remove_duplicate(vector<_outreach>& voutreach);

// automatic load up
bool ProcessPhpLatexCv(void);
// ***************************************************************************
void center_print(uint mode,ostream& oss);

// ***************************************************************************
// references search
void SystemReferences(const string& system_in,vector<string>& vrefs,bool=FALSE);  // if true then put et_al
void SystemReferences(const string& system_in,vector<uint>& vwnumber); 
bool SystemInSystems(const string& system,const string& systems);
bool SystemInSystems(const string& system,const vector<string>& vsystems);
bool AlloyAlphabeticLIBRARY(const string& system);

// ***************************************************************************
// OLD STUFF FOR SECURITY
// bool ProcessSecurityOptions(vector<string> argv,vector<string> cmds);
// void Aflowlib_AUROHOUSE(void);
// void Aflowlib_AUROHOUSE(int deltat);
// void Aflowlib_AUROHOUSE(int deltat,bool CURATOR);
// ***************************************************************************

// ***************************************************************************
//ME20191001
// AflowDB class

namespace aflowlib {

  struct DBStats {
    vector<string> columns;
    vector<vector<int> > count;  // 2D to accommodate bool
    vector<std::pair<string, int> > loop_counts;
    vector<string> max;
    vector<string> min;
    int nentries;
    int nsystems;
    vector<vector<string> > set;
    vector<string> species;
    vector<string> types;
    string catalog;
  };

  class AflowDB : public xStream {
    public:
      AflowDB(const string&, ostream& oss=std::cout);
      AflowDB(const string&, const string&, const string&, ostream& oss=std::cout);
      AflowDB(const AflowDB&);
      AflowDB& operator=(const AflowDB&);
      ~AflowDB();
      void clear();

      bool isTMP();

      int rebuildDatabase(bool force_rebuild=false);
      int rebuildDatabase(const string&, bool force_rebuild=false);
      int rebuildDatabase(const vector<string>&, bool force_rebuild=false);
      int patchDatabase(const string&, bool check_timestamps=false);
      int patchDatabase(const vector<string>&, bool check_timestamps=false);
      void analyzeDatabase(const string&);

      string getEntry(const string&, filetype);
      _aflowlib_entry getEntryAentry(const string&);
      vector<string> getEntrySet(const string&, filetype);
      vector<_aflowlib_entry> getEntrySetAentry(const string&);

      vector<vector<string> > getEntrySetData(const string&);

      vector<string> getTables(const string& where="");
      vector<string> getTables(sqlite3*, const string& where="");

      vector<string> getColumnNames(const string&);
      vector<string> getColumnNames(sqlite3*, const string&);
      vector<string> getColumnTypes(const string&);
      vector<string> getColumnTypes(sqlite3*, const string&);

      vector<vector<string> > getRows(const string&, const string& where="");
      vector<vector<string> > getRows(sqlite3*, const string&, const string& where="");
      vector<vector<string> > getRowsMultiTables(const string& where="");
      vector<vector<string> > getRowsMultiTables(sqlite3*, const string& where="");
      vector<vector<string> > getRowsMultiTables(const vector<string>&, const string& where="");
      vector<vector<string> > getRowsMultiTables(sqlite3*, const vector<string>&, const string& where="");
      string getValue(const string&, const string&, const string& where="");
      string getValue(sqlite3*, const string&, const string&, const string& where="");
      string getProperty(const string&, const string&, const string&, const string& where="");
      string getProperty(sqlite3*, const string&, const string&, const string&, const string& where="");
      vector<string> getPropertyMultiTables(const string&, const vector<string>&, const string&, const string& where="");
      vector<string> getPropertyMultiTables(sqlite3*, const string&, const vector<string>&, const string&, const string& where="");
      vector<string> getSet(const string&, const string&, bool distinct=false, const string& where="", int limit=0, const string& order_by="");
      vector<string> getSet(sqlite3*, const string&, const string&, bool distinct=false, const string& where="", int limit=0, const string& order_by="");
      vector<string> getSetMultiTables(const vector<string>&, const string&, bool distinct=false, const string& where="", int limit=0);
      vector<string> getSetMultiTables(sqlite3*, const vector<string>&, const string&, bool distinct=false, const string& where="", int limit=0);

      void transaction(bool);

    private:
      void free();
      void open(int = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_FULLMUTEX);
      void close();
      void copy(const AflowDB&);
      void initializeExtraSchema();

      sqlite3* db;
      bool is_tmp;
      aurostd::xoption vschema_extra;
      string data_path;
      string database_file;
      string lock_file;

      void openTmpFile(int open_flags=SQLITE_OPEN_READWRITE|SQLITE_OPEN_CREATE|SQLITE_OPEN_FULLMUTEX, bool copy_original=false);
      bool closeTmpFile(bool force_copy=false, bool keep=false, bool nocopy=false);

      void rebuildDB();
      void buildTables(int, int, const vector<string>&, const vector<string>&);
      void populateTable(const string&, const vector<string>&, const vector<string>&, const vector<vector<string> >&);

      uint applyPatchFromJsonl(const vector<string>&);
      bool auidInDatabase(const string&);
      void updateEntry(const string&, const vector<string>&, const vector<string>&);

      vector<string> getSchemaKeys();
      vector<string> getDataTypes(const vector<string>&, bool);
      vector<string> getDataValues(const string&, const vector<string>&, const vector<string>&);

      DBStats initDBStats(const string&, const vector<string>&);
      DBStats getCatalogStats(const string&, const vector<string>&, const vector<string>&);
      void getColStats(int, int, const string&, const vector<string>&, const vector<string>&,
          const vector<string>&, const vector<string>&, vector<vector<vector<int> > >&, vector<vector<int> >&,
          vector<vector<vector<string> > >&, vector<vector<vector<string> > >&);
      vector<string> getUniqueFromJsonArrays(const vector<string>&);
      string stats2json(const DBStats&);

      void createIndex(const string&, const string&, const string&);
      void dropIndex(const string&);
      void dropTable(const string&);
      void createTable(const string&, const vector<string>&, const string&);
      void createTable(const string&, const vector<string>&, const vector<string>&);
      void insertValues(const string&, const vector<string>&);
      void insertValues(const string&, const vector<string>&, const vector<string>&);
      void updateRow(const string&, const vector<string>&, const vector<string>&, const string&);
      string prepareSELECT(const string&, const string&, const string&, const string& where="", int limit=0, const string& order_by="");
      string prepareSELECT(const string&, const string&, const vector<string>&, const string& where="", int limit=0, const string& order_by="");
  };
}  // namespace aflowlib

// ----------------------------------------------------------------------------
// aflowlib_libraries_scrubber.cpp

// will be moved near LI2RAW
namespace aflowlib {
  uint MOSFET(int mode,bool VERBOSE);
  uint MAIL2SCAN(string library,bool VERBOSE);
  uint LIB2SCRUB(string library,bool VERBOSE);
  bool LIB2AUID(string entry,bool TEST,bool _VERBOSE);
}  // namespace aflowlib



#endif //  _AFLOWLIB_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2021           *
// *                                                                         *
// ***************************************************************************
