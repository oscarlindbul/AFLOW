// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2021           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - David Hicks - 2018
// FILE "ANRL/aflow_anrl_AB2C_oC16_40_a_2b_b.cpp"

#ifndef _AFLOW_ANRL_AB2C_oC16_40_a_2b_b_CPP // AFLOW_REMOVE_GREP
#define _AFLOW_ANRL_AB2C_oC16_40_a_2b_b_CPP // AFLOW_REMOVE_GREP
#include "../aflow.h" // AFLOW_REMOVE_GREP

namespace anrl {
  uint WebANRL_AB2C_oC16_40_a_2b_b(stringstream &web,bool LDEBUG);
}

namespace anrl {
  uint PrototypeANRL_AB2C_oC16_40_a_2b_b(stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG) {
    // system AB2C_oC16_40_a_2b_b

    if(XHOST.vflag_control.flag("WWW")) {
      WebANRL_AB2C_oC16_40_a_2b_b(web,LDEBUG); // PLUG WEB STUFF
#ifdef _ANRL_NOWEB_
      web << "no web";
      cout << web.str() << endl;
#else
      cout << web.str() << endl;
#endif
      return 0; //DX20200727
    }

    vector<double> vparameters;
    aurostd::string2tokens(parameters,vparameters,",");

    uint nspecies,natoms,spacegroup,nunderscores,nparameters;
    string label,Pearson_symbol,params,Strukturbericht,prototype,dialect;

    anrl::vproto2tokens(proto_line,label,nspecies,natoms,spacegroup,nunderscores,nparameters,Pearson_symbol,params,Strukturbericht,prototype,dialect);

    anrl::PrototypeANRL_Consistency(vparameters.size(),nparameters,prototype,label,
        Strukturbericht,Pearson_symbol,spacegroup,params,print_mode);    

    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: FOUND" << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: label=" << label << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: nspecies=" << nspecies << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: natoms=" << natoms << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: spacegroup=" << spacegroup << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: nunderscores=" << nunderscores << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: nparameters=" <<  nparameters << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: Pearson_symbol=" << Pearson_symbol << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: params=" << params << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: Strukturbericht=" << Strukturbericht << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: prototype=" << prototype << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: dialect=" << dialect << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: vparameters.size()=" << vparameters.size() << endl;}

    xvector<double> xn(3);   xn(1)=1.0;xn(2)=0.0;xn(3)=0.0;
    xvector<double> yn(3);   yn(1)=0.0;yn(2)=1.0;yn(3)=0.0;
    xvector<double> zn(3);   zn(1)=0.0;zn(2)=0.0;zn(3)=1.0;
    xvector<double> a1(3),a2(3),a3(3);

    if(print_mode==1 && vparameters.size()==0){
      for(uint n=0;n<nparameters;n++){
        vparameters.push_back(0);
      }
    }

    uint i=0;
    double a=vparameters.at(i++);                  if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: a=" << a << endl;}
    double bovera=vparameters.at(i++),b=bovera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: b=" << b << " (b/a=" << bovera << ")" << endl;}
    double covera=vparameters.at(i++),c=covera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: c=" << c << " (c/a=" << covera << ")" << endl;}

    double z1=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: z1=" << z1 << endl;}
    double y2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: y2=" << y2 << endl;}
    double z2=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: z2=" << z2 << endl;}
    double y3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: y3=" << y3 << endl;}
    double z3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: z3=" << z3 << endl;}
    double y4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: y4=" << y4 << endl;}
    double z4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_AB2C_oC16_40_a_2b_b: z4=" << z4 << endl;}

    str.iomode=IOVASP_AUTO;
    str.title=label+" params="+parameters+" SG="+aurostd::utype2string(spacegroup)+DOI_ANRL; //CO20190520
    str.scale=1.0;

    a1=a*xn;
    a2=(1.0/2.0)*b*yn-(1.0/2.0)*c*zn;
    a3=(1.0/2.0)*b*yn+(1.0/2.0)*c*zn;

    str.lattice(1,1)=a1(1);str.lattice(1,2)=a1(2);str.lattice(1,3)=a1(3);
    str.lattice(2,1)=a2(1);str.lattice(2,2)=a2(2);str.lattice(2,3)=a2(3);
    str.lattice(3,1)=a3(1);str.lattice(3,2)=a3(2);str.lattice(3,3)=a3(3);

    // symbolic representation of lattice vectors
    vector<string> a1_equation, a2_equation, a3_equation;
    a1_equation.push_back("a");a1_equation.push_back("0");a1_equation.push_back("0");
    a2_equation.push_back("0");a2_equation.push_back("(1.0/2.0)*b");a2_equation.push_back("-(1.0/2.0)*c");
    a3_equation.push_back("0");a3_equation.push_back("(1.0/2.0)*b");a3_equation.push_back("(1.0/2.0)*c");
    str.symbolic_math_lattice.push_back(a1_equation);
    str.symbolic_math_lattice.push_back(a2_equation);
    str.symbolic_math_lattice.push_back(a3_equation);

    str.num_lattice_parameters = 3;

    str.num_parameters = vparameters.size();
    vector<string> parameter_list; aurostd::string2tokens(params,parameter_list,",");
    str.prototype_parameter_list = parameter_list;
    str.prototype_parameter_values = vparameters;

    if(print_mode!=1){
      str.FixLattices(); // Reciprocal/f2c/c2f
    }

    _atom atom;

    atom.name="A"; atom.type=0;                                       // atom B1
    atom.fpos(1)=0.0;atom.fpos(2)=-z1;atom.fpos(3)=z1;                     // atom B1
    atom.fpos_equation.clear();atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("-z1");atom.fpos_equation.push_back("z1");// atom B1 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B1 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B1

    atom.name="A"; atom.type=0;                                       // atom B2
    atom.fpos(1)=(1.0/2.0);atom.fpos(2)=-z1;atom.fpos(3)=z1;                     // atom B2
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("-z1");atom.fpos_equation.push_back("z1");// atom B2 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B2 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B2

    atom.name="B"; atom.type=1;                                       // atom B3
    atom.fpos(1)=(1.0/4.0);atom.fpos(2)=(y2-z2);atom.fpos(3)=(y2+z2);                     // atom B3
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("(y2-z2)");atom.fpos_equation.push_back("(y2+z2)");// atom B3 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B3 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B3

    atom.name="B"; atom.type=1;                                       // atom B4
    atom.fpos(1)=(3.0/4.0);atom.fpos(2)=(-y2-z2);atom.fpos(3)=(-y2+z2);                     // atom B4
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("(-y2-z2)");atom.fpos_equation.push_back("(-y2+z2)");// atom B4 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B4 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B4

    atom.name="B"; atom.type=1;                                       // atom B5
    atom.fpos(1)=(1.0/4.0);atom.fpos(2)=(y3-z3);atom.fpos(3)=(y3+z3);                     // atom B5
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("(y3-z3)");atom.fpos_equation.push_back("(y3+z3)");// atom B5 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B5 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B5

    atom.name="B"; atom.type=1;                                       // atom B6
    atom.fpos(1)=(3.0/4.0);atom.fpos(2)=(-y3-z3);atom.fpos(3)=(-y3+z3);                     // atom B6
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("(-y3-z3)");atom.fpos_equation.push_back("(-y3+z3)");// atom B6 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B6 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B6

    atom.name="C"; atom.type=2;                                       // atom B7
    atom.fpos(1)=(1.0/4.0);atom.fpos(2)=(y4-z4);atom.fpos(3)=(y4+z4);                     // atom B7
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/4.0)");atom.fpos_equation.push_back("(y4-z4)");atom.fpos_equation.push_back("(y4+z4)");// atom B7 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B7 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B7

    atom.name="C"; atom.type=2;                                       // atom B8
    atom.fpos(1)=(3.0/4.0);atom.fpos(2)=(-y4-z4);atom.fpos(3)=(-y4+z4);                     // atom B8
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(3.0/4.0)");atom.fpos_equation.push_back("(-y4-z4)");atom.fpos_equation.push_back("(-y4+z4)");// atom B8 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B8 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B8


    return str.atoms.size();  
  }
} // namespace anrl

namespace anrl {
  uint WebANRL_AB2C_oC16_40_a_2b_b(stringstream& web,bool LDEBUG) {
#ifndef _ANRL_NOWEB_
#endif

    if(LDEBUG) {cerr << "anrl:: WebANRL_AB2C_oC16_40_a_2b_b: web.str().size()=" << web.str().size() << endl;}

    return web.str().size();
  }
} // namespace anrl

#endif // AFLOW_REMOVE_GREP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2021           *
// *                                                                         *
// ***************************************************************************

