// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2021           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - David Hicks - 2020
// FILE "ANRL/aflow_anrl_A3BC4_oC16_65_ai_b_q.cpp"

#ifndef _AFLOW_ANRL_A3BC4_oC16_65_ai_b_q_CPP // AFLOW_REMOVE_GREP
#define _AFLOW_ANRL_A3BC4_oC16_65_ai_b_q_CPP // AFLOW_REMOVE_GREP
#include "../aflow.h" // AFLOW_REMOVE_GREP

namespace anrl {
  uint WebANRL_A3BC4_oC16_65_ai_b_q(stringstream &web,bool LDEBUG);
}

namespace anrl {
  uint PrototypeANRL_A3BC4_oC16_65_ai_b_q(stringstream &web,xstructure& str,string parameters,string proto_line,uint print_mode,bool LDEBUG) {
    // system A3BC4_oC16_65_ai_b_q

    if(XHOST.vflag_control.flag("WWW")) {
      WebANRL_A3BC4_oC16_65_ai_b_q(web,LDEBUG); // PLUG WEB STUFF
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
        Strukturbericht,Pearson_symbol,spacegroup, params,print_mode);    

    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3BC4_oC16_65_ai_b_q: FOUND" << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3BC4_oC16_65_ai_b_q: label=" << label << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3BC4_oC16_65_ai_b_q: nspecies=" << nspecies << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3BC4_oC16_65_ai_b_q: natoms=" << natoms << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3BC4_oC16_65_ai_b_q: spacegroup=" << spacegroup << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3BC4_oC16_65_ai_b_q: nunderscores=" << nunderscores << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3BC4_oC16_65_ai_b_q: nparameters=" <<  nparameters << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3BC4_oC16_65_ai_b_q: Pearson_symbol=" << Pearson_symbol << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3BC4_oC16_65_ai_b_q: params=" << params << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3BC4_oC16_65_ai_b_q: Strukturbericht=" << Strukturbericht << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3BC4_oC16_65_ai_b_q: prototype=" << prototype << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3BC4_oC16_65_ai_b_q: dialect=" << dialect << endl;}
    if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3BC4_oC16_65_ai_b_q: vparameters.size()=" << vparameters.size() << endl;}

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
    double a=vparameters.at(i++);                  if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3BC4_oC16_65_ai_b_q: a=" << a << endl;}
    double bovera=vparameters.at(i++),b=bovera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3BC4_oC16_65_ai_b_q: b=" << b << " (b/a=" << bovera << ")" << endl;}
    double covera=vparameters.at(i++),c=covera*a;  if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3BC4_oC16_65_ai_b_q: c=" << c << " (c/a=" << covera << ")" << endl;}

    double y3=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3BC4_oC16_65_ai_b_q: y3=" << y3 << endl;}
    double x4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3BC4_oC16_65_ai_b_q: x4=" << x4 << endl;}
    double y4=vparameters.at(i++);                 if(LDEBUG) { cerr << "anrl::PrototypeANRL_A3BC4_oC16_65_ai_b_q: y4=" << y4 << endl;}

    str.iomode=IOVASP_AUTO;
    str.title=label+" params="+parameters+" SG="+aurostd::utype2string(spacegroup)+DOI_ANRL;
    str.scale=1.0;

    a1=(1.0/2.0)*a*xn-(1.0/2.0)*b*yn;
    a2=(1.0/2.0)*a*xn+(1.0/2.0)*b*yn;
    a3=c*zn;

    str.lattice(1,1)=a1(1);str.lattice(1,2)=a1(2);str.lattice(1,3)=a1(3);
    str.lattice(2,1)=a2(1);str.lattice(2,2)=a2(2);str.lattice(2,3)=a2(3);
    str.lattice(3,1)=a3(1);str.lattice(3,2)=a3(2);str.lattice(3,3)=a3(3);

    // symbolic representation of lattice vectors
    vector<string> a1_equation, a2_equation, a3_equation;
    a1_equation.push_back("(1.0/2.0)*a");a1_equation.push_back("-(1.0/2.0)*b");a1_equation.push_back("0");
    a2_equation.push_back("(1.0/2.0)*a");a2_equation.push_back("(1.0/2.0)*b");a2_equation.push_back("0");
    a3_equation.push_back("0");a3_equation.push_back("0");a3_equation.push_back("c");
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
    atom.fpos(1)=0.0;atom.fpos(2)=0.0;atom.fpos(3)=0.0;                     // atom B1
    atom.fpos_equation.clear();atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("0.0");atom.fpos_equation.push_back("0.0");// atom B1 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B1 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B1

    atom.name="A"; atom.type=0;                                       // atom B3
    atom.fpos(1)=-y3;atom.fpos(2)=y3;atom.fpos(3)=0.0;                     // atom B3
    atom.fpos_equation.clear();atom.fpos_equation.push_back("-y3");atom.fpos_equation.push_back("y3");atom.fpos_equation.push_back("0.0");// atom B3 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B3 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B3

    atom.name="A"; atom.type=0;                                       // atom B4
    atom.fpos(1)=y3;atom.fpos(2)=-y3;atom.fpos(3)=0.0;                     // atom B4
    atom.fpos_equation.clear();atom.fpos_equation.push_back("y3");atom.fpos_equation.push_back("-y3");atom.fpos_equation.push_back("0.0");// atom B4 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B4 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B4

    atom.name="B"; atom.type=1;                                       // atom B2
    atom.fpos(1)=(1.0/2.0);atom.fpos(2)=(1.0/2.0);atom.fpos(3)=0.0;                     // atom B2
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("(1.0/2.0)");atom.fpos_equation.push_back("0.0");// atom B2 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B2 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B2

    atom.name="C"; atom.type=2;                                       // atom B5
    atom.fpos(1)=(x4-y4);atom.fpos(2)=(x4+y4);atom.fpos(3)=(1.0/2.0);                     // atom B5
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x4-y4)");atom.fpos_equation.push_back("(x4+y4)");atom.fpos_equation.push_back("(1.0/2.0)");// atom B5 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B5 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B5

    atom.name="C"; atom.type=2;                                       // atom B6
    atom.fpos(1)=(-x4+y4);atom.fpos(2)=(-x4-y4);atom.fpos(3)=(1.0/2.0);                     // atom B6
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x4+y4)");atom.fpos_equation.push_back("(-x4-y4)");atom.fpos_equation.push_back("(1.0/2.0)");// atom B6 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B6 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B6

    atom.name="C"; atom.type=2;                                       // atom B7
    atom.fpos(1)=(-x4-y4);atom.fpos(2)=(-x4+y4);atom.fpos(3)=(1.0/2.0);                     // atom B7
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(-x4-y4)");atom.fpos_equation.push_back("(-x4+y4)");atom.fpos_equation.push_back("(1.0/2.0)");// atom B7 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B7 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B7

    atom.name="C"; atom.type=2;                                       // atom B8
    atom.fpos(1)=(x4+y4);atom.fpos(2)=(x4-y4);atom.fpos(3)=(1.0/2.0);                     // atom B8
    atom.fpos_equation.clear();atom.fpos_equation.push_back("(x4+y4)");atom.fpos_equation.push_back("(x4-y4)");atom.fpos_equation.push_back("(1.0/2.0)");// atom B8 // symbolic math for atom positions
    str.comp_each_type.at(atom.type)+=1.0;                            // atom B8 // if we need partial occupation
    str.atoms.push_back(atom);                                        // atom B8


    return str.atoms.size();  
  }
} // namespace anrl

namespace anrl {
  uint WebANRL_A3BC4_oC16_65_ai_b_q(stringstream& web,bool LDEBUG) {
#ifndef _ANRL_NOWEB_
#endif

    if(LDEBUG) {cerr << "anrl:: WebANRL_A3BC4_oC16_65_ai_b_q: web.str().size()=" << web.str().size() << endl;}

    return web.str().size();
  }
} // namespace anrl

#endif // AFLOW_REMOVE_GREP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2021           *
// *                                                                         *
// ***************************************************************************

