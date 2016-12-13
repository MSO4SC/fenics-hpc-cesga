#ifndef NSESOLVER_H
#define NSESOLVER_H

#include <dolfin.h>
#include <dolfin/main/MPI.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/fem/UFC.h>
#ifdef ENABLE_UFL 
#include "ufc2/NSEMomentum3D.h"
#include "ufc2/NSEContinuity3D.h"
#include "ufc2/NSEDualMomentum3D.h"
#include "ufc2/NSEDualContinuity3D.h"
#include "ufc2/NSEErrRepMomentum3D.h"
#include "ufc2/NSEErrRepContinuity3D.h"
#include "ufc2/Drag3D.h"
#include "ufc2/NSEH1.h"
#include "ufc2/NSEH12.h"
#include "ufc2/NSEH1Momentum3D.h"
#include "ufc2/NSEH1Continuity3D.h"
#include "ufc2/NSEH1MomentumGlobal3D.h"
#include "ufc2/NSEH1ContinuityGlobal3D.h"
#include "ufc2/NSEMomentumResidual3D.h"
#include "ufc2/NSEContinuityResidual3D.h"
#include "ufc2/NSEMomentumResidualGlobal3D.h"
#include "ufc2/NSEContinuityResidualGlobal3D.h"
#include "ufc2/NSEErrEst.h"
#include "ufc2/NSEErrEstGlobal.h"
#else
#include "ufc2/NSEMomentum3D.h"
#include "ufc2/NSEContinuity3D.h"
#endif

#include "dolfin/NodeNormal.h"
#include "dolfin/SpaceTimeFunction.h"
#include "dolfin/SlipBC.h"

#include <ostream>
#include <iomanip>
#include <cstring>
#include <sstream>
#include <string>
#include <algorithm>
#include <map>
#include <mpi.h>

using namespace dolfin;

typedef std::vector<std::pair<SubDomain*,Function*>> DirichletBCList;
typedef std::vector<SubDomain*> SlipBCList;

// Comparison operator for index/value pairs
struct less_pair : public std::binary_function<std::pair<int, real>,
                                               std::pair<int, real>, bool>
{
  bool operator()(std::pair<int, real> x, std::pair<int, real> y)
  {
    return x.second < y.second;
  }
};


class NSESolver
{
public:
    enum class SolverType {primalSolver, dualSolver};

    NSESolver(
        Mesh& mesh,
        DirichletBCList dbcs_m,
        SlipBCList sbcs_m,
        DirichletBCList dbcs_c,
        DirichletBCList dbcs_dm,
        Function* thetadrag,
        Function* thetalift,
        Function* sm,
        Function* psim,
        Function* psic,
        Function* bpsim
        ) :
      mesh(mesh),
      thetadrag(thetadrag),
      thetalift(thetalift),
      sm(sm),
      psim(psim),
      psic(psic),
      bpsim(bpsim)
    {
      real theta = 0.;

//      // Parse command-line arguments
//      if(argc >= 2)
//      {
//        simcase = std::string(argv[1]);
//
//      }
//
//      if(argc >= 3)
//      {
//        theta = atof(argv[2]);
//      }

      dolfin_set("output destination","silent");
      if(dolfin::MPI::processNumber() == 0)
        dolfin_set("output destination","terminal");

      assembler = new Assembler(mesh);

      message("Running on %d %s", dolfin::MPI::numProcesses(), 
              (dolfin::MPI::numProcesses() > 1 ? "nodes" : "node"));
      message("Global number of vertices: %d", 
              (dolfin::MPI::numProcesses() > 1 ? mesh.distdata().global_numVertices() : mesh.numVertices()));
      message("Global number of cells: %d", 
              (dolfin::MPI::numProcesses() > 1 ? mesh.distdata().global_numCells() : mesh.numCells()));

    //    for(int i = 0; i < 1; i++)
    //      mesh.refine();

      message("Global number of vertices: %d", 
              (dolfin::MPI::numProcesses() > 1 ? mesh.distdata().global_numVertices() : mesh.numVertices()));
      message("Global number of cells: %d", 
              (dolfin::MPI::numProcesses() > 1 ? mesh.distdata().global_numCells() : mesh.numCells()));

      h = new MeshSize(mesh);
#warning "unused variable"
      n = new FacetNormal(mesh);
      nn = new NodeNormal(mesh);
      cv = new CellVolume(mesh);

      // create bcs
      for (auto bc : dbcs_m)
      {
        bcs_m.push_back(new DirichletBC(*bc.second,mesh,*bc.first));
      }

      for (auto bc : sbcs_m)
      {
        bcs_m.push_back(new SlipBC(mesh,*bc,*nn));
      }

      for (auto bc : dbcs_dm)
      {
        bcs_dm.push_back(new DirichletBC(*bc.second,mesh,*bc.first));
      }

#warning "this is a temporary ugly thing, please improve me!"
      dbc_c = new DirichletBC(*dbcs_c[0].second,mesh,*dbcs_c[0].first);

#warning "unused variables"
      uint *c_indices = 0;
      uint *indices = 0;

      hmin = h->min();
      cout << "hmin: " << hmin << endl;
      
      k = 0.1*hmin;
      nu = 0.0;
      c1 = 0.1;
      c2 = 0.0;
      c3 = 0.0;

      rdtol = 5e-2;
      maxit = 10;

      cout << "c1: " << c1 << endl;
      cout << "c2: " << c2 << endl;

      // Declare all needed functions
      u = new Function;
      u0 = new Function;
      p = new Function;
      p0 = new Function;
      nuf = new Function(mesh, nu);
      kf = new Function(mesh, k);
      c1f = new Function(mesh, c1);
      c2f = new Function(mesh, c2);
      hminf = new Function(mesh, c3*hmin);
      dtu = new Function;
#warning "unused variableS"
      Function k0f(mesh, 0.0);
      Function umean;

      Function X;

      rd_u = new Function;
      rd_p = new Function;

      up = new Function;
      pp = new Function;
      up0 = new Function;
      dtup = new Function;

      Rm = new Function;
      Rmtot = new Function;
      Rc = new Function;
      Rctot = new Function;
      wm = new Function;
      wmtot = new Function;
      wc = new Function;
      wctot = new Function;

      ei_m = new Function;
      ei_c = new Function;
      eij_m = new Function;
      eij_c = new Function;
      eif = new Function;

      tau_1 = new Function;
      tau_2 = new Function;
      normal = new Function;

      ei = new Vector;

      ap_m = new NSEMomentum3DBilinearForm(*u, *p, *nuf, *h, *kf, *c1f, *u0,
          *normal, *sm);
      Lp_m = new NSEMomentum3DLinearForm(*u, *p, *nuf, *h, *kf, *c1f, *u0,
          *normal, *sm);

      ap_c = new NSEContinuity3DBilinearForm(*h, *kf, *c1f, *hminf);
      Lp_c = new NSEContinuity3DLinearForm(*u, *p, *h, *kf, *c1f, *u0, *p0,
          *hminf);

      //NSEDualMomentum3DBilinearForm ad_m(*u, *up, *nuf, *h, *kf, *c1f, c2f, *u0);
      ad_m = new NSEDualMomentum3DBilinearForm(*u, *up, *p, *nuf, *h, *kf,
          *c1f, *c2f, *u0);
      Ld_m = new NSEDualMomentum3DLinearForm(*u, *up, *p, *nuf, *h, *kf,
          *c1f, *c2f, *u0, *psim, *bpsim);

      ad_c = new NSEDualContinuity3DBilinearForm(*u, *h, *kf, *c1f, *c2f,
          *hminf);
      //NSEDualContinuity3DLinearForm Ld_c(*u, *up, *p, *h, *kf, *c1f, *u0, *p0, *hminf, psic);
      Ld_c = new NSEDualContinuity3DLinearForm(*u, *p, *h, *kf, *c1f, *c2f,
          *u0, *p0, *hminf, *psic);

      Lrep_m = new NSEErrRepMomentum3DLinearForm(*up, *pp, *nuf, *dtup, *u);
      Lrep_c = new NSEErrRepContinuity3DLinearForm(*up, *p);

      Md = new Drag3DFunctional(*normal, *thetadrag, *p);
      Ml = new Drag3DFunctional(*normal, *thetalift, *p);

      MH1 = new NSEH1Functional(*u, *p, *h);
      MH12 = new NSEH12Functional(*u, *p);
      MHerrest = new NSEErrEstFunctional(*h, *cv, *Rm, *Rc, *wm, *wc);
      MHerrestg = new NSEErrEstGlobalFunctional(*h, *cv, *Rm, *Rc, *wm, *wc);

      LRm = new NSEMomentumResidual3DLinearForm(*u, *p, *kf, *u0);
      MRm = new NSEMomentumResidual3DFunctional(*h, *cv, *Rm);
      LRc = new NSEContinuityResidual3DLinearForm(*u, *u0);
      MRc = new NSEContinuityResidual3DFunctional(*h, *cv, *Rc);
      MRgm = new NSEMomentumResidualGlobal3DFunctional(*u, *p, *h, *kf, *u0);
      MRgc = new NSEContinuityResidualGlobal3DFunctional(*u, *h, *u0);

      Lwm = new NSEH1Momentum3DLinearForm(*u, *kf, *u0);
      Mwm = new NSEH1Momentum3DFunctional(*h, *cv, *wm);
      Lwc = new NSEH1Continuity3DLinearForm(*p);
      Mwc = new NSEH1Continuity3DFunctional(*h, *cv, *wc);

      Mgwm = new NSEH1MomentumGlobal3DFunctional(*u, *h, *kf, *u0);
      Mgwc = new NSEH1ContinuityGlobal3DFunctional(*p, *h);

      // Initialize functions
      u->init(mesh, *ap_m, 0);
      u0->init(mesh, *ap_m, 0);
      p->init(mesh, *ap_c, 0);
      p0->init(mesh, *ap_c, 0);

      X.init(mesh, *ap_m, 0);

      rd_u->init(mesh, *ap_m, 0);
      rd_p->init(mesh, *ap_c, 0);

      up->init(mesh, *ap_m, 0);
      up0->init(mesh, *ap_m, 0);
      pp->init(mesh, *ap_c, 0);
      dtup->init(mesh, *ap_m, 0);

      dtu->init(mesh, *ap_m, 0);

      Rm->init(mesh, *LRm, 0);
      Rc->init(mesh, *LRc, 0);
      wm->init(mesh, *Lwm, 0);
      wc->init(mesh, *Lwc, 0);

      Rmtot->init(mesh, *LRm, 0);
      Rctot->init(mesh, *LRc, 0);
      wmtot->init(mesh, *Lwm, 0);
      wctot->init(mesh, *Lwc, 0);

      ei_m->init(mesh, *Lrep_m, 0);
      eij_m->init(mesh, *Lrep_m, 0);
      ei_c->init(mesh, *Lrep_c, 0);
      eij_c->init(mesh, *Lrep_c, 0);
      eif->init(mesh, *ei, *Lrep_c, 0);

      normal->init(mesh, *ap_m, 0);
      tau_1->init(mesh, *ap_m, 0);
      tau_2->init(mesh, *ap_m, 0);

      p->vector() = 1.0;
      p0->vector() = 1.0;

      ei_m->vector() = 0.0;
      ei_c->vector() = 0.0;

      ComputeTangentialVectors((Vector&)tau_1->vector(),
          (Vector&)tau_2->vector(),
          (Vector&)normal->vector(), *ap_m, *nn);

      dolfin_set("PDE linear solver", "iterative");

      // Declare PDE solvers
      pdep_m = new LinearPDE(*ap_m, *Lp_m, mesh, bcs_m);
      pdep_c = new LinearPDE(*ap_c, *Lp_c, mesh, dbc_c[0], cg);

      pded_m = new LinearPDE(*ad_m, *Ld_m, mesh, bcs_dm);
      pded_c = new LinearPDE(*ad_c, *Ld_c, mesh, dbc_c[0]);

      U = new Function;
      P = new Function;

      computeX(X, ap_m, mesh);
      rotate(X, ap_m, mesh, theta);

      file_u = new File("velocity.bin");
      file_p = new File("pressure.bin");
      file_du = new File("dvelocity.bin");
      file_dp = new File("dpressure.bin");
      file_m = new File("mesh.bin");

      //file_m << mesh;

      File meshfile("mesh_out.bin");
      meshfile << mesh;

      Up = new SpaceTimeFunction(mesh, *up);
      Up->setBasename("velocity_v");

      dtUp = new SpaceTimeFunction(mesh, *dtup);
      dtUp->setBasename("dtvelocity_v");

      Pp = new SpaceTimeFunction(mesh, *pp);
      Pp->setBasename("pressure_v");

      Rmp = new SpaceTimeFunction(mesh, *Rm);
      Rmp->setBasename("Rm_v");

      Rcp = new SpaceTimeFunction(mesh, *Rc);
      Rcp->setBasename("Rc_v");
    }

    virtual ~NSESolver ();

    template<SolverType st>
    void step()
    {
        stimer = time();

        s = primal_T - t;

        umax = u->vector().norm(linf);
        
        step_preSolve<st>();
        
        if(stepcounter >= 100)
        {
          if(iteration0 > 5)
            stabcounter = 10;
          k = 4.0*hmin/std::max(1., umax);
          if(hmin >= 0.1)
          {
            k = 0.5*hmin/std::max(1., umax);
          }
          if(stabcounter > 0)
            k /= 4.;
            kf->init(mesh, k);
        }

        solve_system();

        step_postSolve<st>();

        u0->vector() = u->vector();
        up0->vector() = up->vector();

        t += k;
        stepcounter++;
        if(stabcounter > 0)
          stabcounter--;
    }

    template <SolverType st = SolverType::primalSolver>
    void run_solver()
    {
        run_preStepping<st>();

        u->vector() = 0.0;
        u0->vector() = 0.0;
        p->vector() = 0.0;
        p0->vector() = 0.0;
        up0->vector() = 0.0;

        Rmtot->vector() = 0.0;
        Rctot->vector() = 0.0;
        wmtot->vector() = 0.0;
        wctot->vector() = 0.0;
        
        t = 0;

        k = 0.1*hmin;
        kf->init(mesh, k);

        // Time-stepping
        while(t <= T)
        {
            step<st>();
        }
        
        run_postStepping<st>();

        cout << "Solver done" << endl;
    }

    void run()
    {
      run_solver<SolverType::primalSolver>();
      run_solver<SolverType::dualSolver>();
    }

    real getT() const {return T;}
    void setT(const real _T)
    {
        T = _T;
        primal_T = T;
        dual_T = 1. * T / 2;
    }


private:
    template <SolverType st>
    inline void step_preSolve()
    {};

    template <SolverType st>
    void run_preStepping()
    {};

    template <SolverType st>
    void run_postStepping()
    {};

    template <SolverType st>
    void step_postSolve()
    {};

    void solve_system()
    {
      // Fixed-point iteration
      for(int i = 0; i < maxit; i++)
      {
        cout << "Solving momentum" << endl;
        real timer = time();
        pde_m->solve(*U);
        rd_u->vector() = U->vector();
        rd_u->vector() -= u->vector();
        real rd_u_norm = rd_u->vector().norm(l2);
        u->vector() = U->vector();
        
        cout << "Solving continuity" << endl;
        pde_c->solve(*P);
        p->vector() = P->vector();
        rd_p->vector() = p->vector();
        rd_p->vector() -= p0->vector();
        real rd_p_norm = rd_p->vector().norm(l2);
        
        cout << "Iteration info: " <<
          "Unorm: " << U->vector().norm(linf) <<
          " Pnorm: " << P->vector().norm(linf) <<
          " Uincr: " <<  rd_u_norm <<
          " Pincr: " <<  rd_p_norm <<
          " k: " << k <<
          " step: " << stepcounter <<
          " t: " << t <<
          " timer: " << time() - timer << endl;
        cout << "iteration: " << i << endl;
        iteration0 = i;
        if(rd_u_norm / u->vector().norm(l2) <= rdtol &&
           rd_p_norm / p->vector().norm(l2) <= rdtol)
        {
          cout << "Step info: " <<
            "Unorm: " << U->vector().norm(linf) <<
            " Pnorm: " << P->vector().norm(linf) <<
            " Uincr: " <<  rd_u_norm / u->vector().norm(l2) <<
            " Pincr: " <<  rd_p_norm / p->vector().norm(l2) <<
            " k: " << k <<
            " step: " << stepcounter <<
            " iters: " << iteration0 + 1 <<
            " t: " << t <<
            " timer: " << time() - stimer << endl;
          break;
        }
        p0->vector() = p->vector();
      }
    }

    void ComputeTangentialVectors(Vector& tau_1,
            Vector& tau_2, Vector& normal, Form& form, NodeNormal& node_normal);
    void merge(real *a,real *b,real *res,int an,int bn);
    void ComputeLargestIndicators_cell(Mesh& mesh, Vector& e_indx,
            std::vector<int>& cells, real percentage);
    void ComputeRefinementMarkers(Mesh& mesh, real percentage,
            Vector& e_indx, MeshFunction<bool>& cell_refinement_marker);
    void computeX(Function& XX, Form* aM, Mesh& mesh);
    void rotate(Function& XX, Form* aM, Mesh& mesh, double theta);


    Mesh& mesh;

//    SolverType solver = SolverType::primalSolver;

    bool coeffchanged = true;
    real int_errest_gstcs = 0;

    real k;
    real nu;
    real c1;
    real c2;
    real c3;

    real rdtol = 5e-2;
    int maxit = 10;

    // Declare all needed functions
    Function* u;
    Function* u0;
    Function* p;
    Function* p0;
    Function* nuf;
    Function* kf;
    Function* c1f;
    Function* c2f;
    Function* hminf;

    Function* up;
    Function* pp;
    Function* up0;
    Function* dtu;
    Function* dtup;

    Function *Rm, *Rmtot;
    Function *Rc, *Rctot;
    Function *wm, *wmtot;
    Function *wc, *wctot;  

    // Declare primal and dual forms
    Form *a_m, *L_m, *a_c, *L_c;

    NSEMomentum3DBilinearForm* ap_m;
    NSEMomentum3DLinearForm* Lp_m;

    NSEContinuity3DBilinearForm* ap_c;
    NSEContinuity3DLinearForm* Lp_c;

    NSEDualMomentum3DBilinearForm* ad_m;
    NSEDualMomentum3DLinearForm* Ld_m;

    NSEDualContinuity3DBilinearForm* ad_c;
    NSEDualContinuity3DLinearForm* Ld_c;


    // Declare PDE solvers
    LinearPDE *pde_m, *pde_c;

    LinearPDE* pdep_m;
    LinearPDE* pdep_c;

    LinearPDE* pded_m;
    LinearPDE* pded_c;

    real T = 0;
    real primal_T = 0;
    real dual_T = 0;

    real hmin;

    SpaceTimeFunction* Up = nullptr;
    SpaceTimeFunction* dtUp = nullptr;
    SpaceTimeFunction* Pp = nullptr;
    SpaceTimeFunction* Rmp = nullptr;
    SpaceTimeFunction* Rcp = nullptr;

    int iteration0 = 0;
    int stabcounter = 0;

    int no_samples = 200;

    Function* U;
    Function* P;

    Function* rd_u;
    Function* rd_p;

    Assembler* assembler;

    Function *ei_m;
    Function *ei_c;
    Function *eij_m;
    Function *eij_c;
    Function *eif;

    NSEErrRepMomentum3DLinearForm *Lrep_m;
    NSEErrRepContinuity3DLinearForm *Lrep_c;

    Drag3DFunctional *Md;
    Drag3DFunctional *Ml;

    NSEH1Functional *MH1;
    NSEH12Functional *MH12;
    NSEErrEstFunctional *MHerrest;
    NSEErrEstGlobalFunctional *MHerrestg;

    NSEMomentumResidual3DLinearForm *LRm;
    NSEMomentumResidual3DFunctional *MRm;
    NSEContinuityResidual3DLinearForm *LRc;
    NSEContinuityResidual3DFunctional *MRc;
    NSEMomentumResidualGlobal3DFunctional *MRgm;
    NSEContinuityResidualGlobal3DFunctional *MRgc;

    NSEH1Momentum3DLinearForm *Lwm;
    NSEH1Momentum3DFunctional *Mwm;
    NSEH1Continuity3DLinearForm *Lwc;
    NSEH1Continuity3DFunctional *Mwc;

    NSEH1MomentumGlobal3DFunctional *Mgwm;
    NSEH1ContinuityGlobal3DFunctional *Mgwc;

    File *file_u;
    File *file_p;
    File *file_du;
    File *file_dp;
    File *file_m;

    MeshSize *h;
    FacetNormal *n;
    NodeNormal *nn;
    CellVolume* cv;
    Vector* ei;

    Function *tau_1;
    Function *tau_2;
    Function *normal;

    // Create boundary conditions
    Function *thetadrag;
    Function *thetalift;
    Function *sm;

    DirichletBC* dbc_c;
    Array<BoundaryCondition*> bcs_m;
    Array<BoundaryCondition*> bcs_dm;

    Function *psim;
    Function *psic;
    Function *bpsim;

    real adapt_percent = 5.;

    int stepcounter = 0;
    real t = 0;
    real s = 0;
    real stimer;
    real umax;

    real mean_drag = 0;
    real tot_lift = 0;
    int n_mean = 0;

    real tot_Rgstm = 0;
    real tot_Rgstc = 0;
    real tot_H1primal = 0;
    real tot_H1primal2 = 0;
    real tot_H1dualm = 0;
    real tot_H1dualc = 0;
    real tot_H1dualgm = 0;
    real tot_H1dualgc = 0;
    real tot_H1dualgstm = 0;
    real tot_H1dualgstc = 0;
    real tot_Rm = 0;
    real tot_Rc = 0;
    real tot_Rgm = 0;
    real tot_Rgc = 0;

    real int_errest_cs = 0;
    real int_errest_gcs = 0;
    
    int sample = 0;
};


template<>
void NSESolver::run_preStepping<NSESolver::SolverType::primalSolver>()
{
  cout << "Starting primal solver" << endl;
  a_m = ap_m; L_m = Lp_m; a_c = ap_c; L_c = Lp_c;
  pde_m = pdep_m; pde_c = pdep_c;
}

template<>
void NSESolver::run_preStepping<NSESolver::SolverType::dualSolver>()
{
  c1 = 4.0;
  c1f->init(mesh, c1);

  cout << "Starting dual solver" << endl;
  a_m = ad_m; L_m = Ld_m; a_c = ad_c; L_c = Ld_c;
  pde_m = pded_m; pde_c = pded_c;
  T = dual_T;
}

template <>
inline void NSESolver::step_preSolve<NSESolver::SolverType::dualSolver>()
{
    cout << "eval dual" << endl;
    Up->eval(s);
    dtUp->eval(s);
    Pp->eval(s);
    Rmp->eval(s);
    Rcp->eval(s);
    cout << "eval dual done" << endl;

    umax = up->vector().norm(linf);
}

template <>
void NSESolver::run_postStepping<NSESolver::SolverType::primalSolver>()
{
#warning "poorly named variable"
  cout << "mean drag: " << mean_drag << endl;
  cout << "total H1primal: " << sqrt(tot_H1primal) << endl;
  cout << "total H1primal2: " << sqrt(tot_H1primal2) << endl;
  cout << "total Rgstm: " << sqrt(tot_Rgstm) << endl;
  cout << "total Rgstc: " << sqrt(tot_Rgstc) << endl;

  int_errest_gstcs = (sqrt(tot_Rgstm) + sqrt(tot_Rgstc));
}

template <>
void NSESolver::run_postStepping<NSESolver::SolverType::dualSolver>()
{
      cout << "Preparing adaptivity" << endl;
      // Adaptive error control
      if(!ParameterSystem::parameters.defined("adapt_algorithm"))
        dolfin_add("adapt_algorithm", "rivara");
      dolfin_set("adapt_algorithm", "rivara");
      if(!ParameterSystem::parameters.defined("output_format"))
        dolfin_add("output_format", "binary");
      dolfin_set("output_format", "binary");
      MeshFunction<bool> cell_marker;

      
      eif->vector() = 0.0;
      eif->vector() += ei_m->vector();
      eif->vector() += ei_c->vector();
      //eif->vector() /= dual_T;
      
      NSEErrRepMomentum3DFunctional M_ei(*eif, *cv);
      real errest = fabs(assembler->assemble(M_ei));

      //int_errest_cs = sqrt(int_errest_cs);

      int_errest_gstcs *= (sqrt(tot_H1dualgstm) + sqrt(tot_H1dualgstc));

      cout << "error estimate: " << errest / dual_T << endl;
      cout << "error estimate cs: " << int_errest_cs / dual_T << endl;
      cout << "error estimate gcs: " << int_errest_gcs / dual_T << endl;
      cout << "error estimate gstcs: " << int_errest_gstcs / dual_T << endl;
      cout << "total H1dualm: " << tot_H1dualm << endl;
      cout << "total H1dualc: " << tot_H1dualc << endl;
      cout << "total H1dualgm: " << tot_H1dualgm << endl;
      cout << "total H1dualgc: " << tot_H1dualgc << endl;
      cout << "total H1dualgstm: " << sqrt(tot_H1dualgstm) << endl;
      cout << "total H1dualgstc: " << sqrt(tot_H1dualgstc) << endl;
      cout << "total Rm: " << tot_Rm << endl;
      cout << "total Rc: " << tot_Rc << endl;
      cout << "total Rgm: " << tot_Rgm << endl;
      cout << "total Rgc: " << tot_Rgc << endl;

      File file_ei("ei.bin");
      file_ei << *eif;

      File file_Rmtot("Rmtot.bin");
      file_Rmtot << *Rmtot;
      File file_Rctot("Rctot.bin");
      file_Rctot << *Rctot;
      File file_wmtot("wmtot.bin");
      file_wmtot << *wmtot;
      File file_wctot("wctot.bin");
      file_wctot << *wctot;

      MeshFunction<real> eimf;
      eimf.init(mesh, mesh.topology().dim());
  
      // Initialize eimf - assumption on dofmap for DG0
      for (CellIterator c(mesh); !c.end(); ++c)
      {
        eimf.set(*c, eif->vector()[c->index()]);
      }

      MPI_Barrier(dolfin::MPI::DOLFIN_COMM);
  
      cout << "Output eimf: " << endl;
      File file_eimf("eimf.bin");
      file_eimf << eimf;
      

      ComputeRefinementMarkers(mesh, adapt_percent, *ei, cell_marker);

      if(MPI::processNumber() == 0)
        dolfin_set("output destination","terminal");
      message("Adaptive refinement");
      message("cells before: %d",
              (dolfin::MPI::numProcesses() > 1 ?
               mesh.distdata().global_numCells() : mesh.numCells()));
      message("vertices before: %d",
              (dolfin::MPI::numProcesses() > 1 ?
               mesh.distdata().global_numVertices() : mesh.numVertices()));
      dolfin_set("output destination","silent");

      RivaraRefinement::refine(mesh, cell_marker);

      if(MPI::processNumber() == 0)
        dolfin_set("output destination","terminal");
      message("cells after: %d",
              (dolfin::MPI::numProcesses() > 1 ?
               mesh.distdata().global_numCells() : mesh.numCells()));
      message("vertices after: %d",
              (dolfin::MPI::numProcesses() > 1 ?
               mesh.distdata().global_numVertices() : mesh.numVertices()));
      dolfin_set("output destination","silent");
      
      File file_rm("rmesh.bin");
      file_rm << mesh;
}

template <>
void NSESolver::step_postSolve<NSESolver::SolverType::primalSolver>()
{
  real drag = 0.0, lift = 0.0;
  drag = assembler->assemble(*Md);
  cout << "drag: " << drag << " t = " << t << endl;
  lift = assembler->assemble(*Ml);
  cout << "lift: " << lift << " t = " << t << endl;

  assembler->assemble(Rm->vector(), *LRm);
  real Rmi = assembler->assemble(*MRm);
  cout << "step primal Rm: " << Rmi << endl;
  assembler->assemble(Rc->vector(), *LRc);
  real Rci = assembler->assemble(*MRc);
  cout << "step primal Rc: " << Rci << endl;

  if(t >= dual_T)
  {
    // Output drag and lift, together with other diagnostics
    mean_drag = (drag + n_mean*mean_drag) / (n_mean + 1);
    cout << "step t: " << t <<
      " drag: " << drag <<
      " lift: " << lift << endl;
    real H1primal = assembler->assemble(*MH1);
    tot_H1primal = (H1primal + n_mean*tot_H1primal) / (n_mean + 1);
    real H1primal2 = assembler->assemble(*MH12);
    tot_H1primal2 = (H1primal2 + n_mean*tot_H1primal2) / (n_mean + 1);
    cout << "step H1 primal: " << tot_H1primal << endl;
    cout << "step H1 primal2: " << tot_H1primal2 << endl;
    n_mean++;

    real Rgstmi = assembler->assemble(*MRgm);
    tot_Rgstm += k*Rgstmi;
    real Rgstci = assembler->assemble(*MRgc);
    tot_Rgstc += k*Rgstci;
  }

  if(stepcounter == 0 || t > T*(real(sample)/real(no_samples)))
  {
    *file_u << *u;
    *file_p << *p;

    // Save primal velocity
    up->vector() = u->vector(); 
    up->vector() += u0->vector(); 
    up->vector() /= 2.;
    File ubinfile(Up->getNewFilename(t));
    ubinfile << up->vector();

    // Save primal velocity time-derivative
    dtu->vector() = u->vector();
    dtu->vector() -= u0->vector();
    dtu->vector() /= k;
    File dtubinfile(dtUp->getNewFilename(t));
    dtubinfile << dtu->vector();

    // Save primal pressure
    File pbinfile(Pp->getNewFilename(t));
    pbinfile << p->vector();

    // Save primal residuals
    File Rmbinfile(Rmp->getNewFilename(t));
    Rmbinfile << Rm->vector();

    File Rcbinfile(Rcp->getNewFilename(t));
    Rcbinfile << Rc->vector();
  }
  else
  {
    *file_du << *u; *file_dp << *p;
  }
  
  sample++;
}

template <>
void NSESolver::step_postSolve<NSESolver::SolverType::dualSolver>()
{
    cout << "errest" << endl;
    assembler->assemble(eij_m->vector(), *Lrep_m);
    ei_m->vector().axpy(k, eij_m->vector());
    assembler->assemble(eij_c->vector(), *Lrep_c);
    ei_c->vector().axpy(k, eij_c->vector());
    cout << "errest done: " << ei_m->vector().norm(linf) <<
      " " << ei_c->vector().norm(linf) << endl;
    
    assembler->assemble(wm->vector(), *Lwm);
    real H1dualm = assembler->assemble(*Mwm);
    tot_H1dualm += k*H1dualm;
    real H1dualgm = sqrt(assembler->assemble(*Mgwm));
    tot_H1dualgm += k*H1dualgm;
    real H1dualgstm = assembler->assemble(*Mgwm);
    tot_H1dualgstm += k*H1dualgstm;
    wmtot->vector().axpy(k, wm->vector());
    assembler->assemble(wc->vector(), *Lwc);
    real H1dualc = assembler->assemble(*Mwc);
    tot_H1dualc += k*H1dualc;
    real H1dualgc = sqrt(assembler->assemble(*Mgwc));
    tot_H1dualgc += k*H1dualgc;
    real H1dualgstc = assembler->assemble(*Mgwc);
    tot_H1dualgstc += k*H1dualgstc;
    wctot->vector().axpy(k, wc->vector());
    
    real Rmi = assembler->assemble(*MRm);
    tot_Rm += k*Rmi;
    real Rgmi = 0.0;
    Rmtot->vector().axpy(k, Rm->vector());
    real Rci = assembler->assemble(*MRc);
    tot_Rc += k*Rci;
    real Rgci = 0.0;
    Rctot->vector().axpy(k, Rc->vector());
    
    real errest_cs = assembler->assemble(*MHerrest);
    int_errest_cs += k*errest_cs;
    real errest_gcs = sqrt(assembler->assemble(*MHerrestg));
    int_errest_gcs += k*errest_gcs;
    real errest_gstcs = assembler->assemble(*MHerrestg);
    int_errest_gstcs += k*errest_gstcs;
    n_mean++;
    
    cout << "step dual t: " << t <<
      " dualm: " << H1dualm <<
      " dualc: " << H1dualc <<
      " dualgm: " << H1dualgm <<
      " dualgc: " << H1dualgc <<
      " Rm: " << Rmi <<
      " Rc: " << Rci <<
      " Rgm: " << Rgmi <<
      " Rgc: " << Rgci <<
      " errest_cs: " << errest_cs <<
      " errest_gcs: " << errest_gcs <<
      endl;
}

void NSESolver::ComputeTangentialVectors(Vector& tau_1, 
                              Vector& tau_2, Vector& normal,
                              Form& form, NodeNormal& node_normal)
{
  UFC ufc(form.form(), mesh, form.dofMaps());
  Cell c(mesh, 0);
  uint local_dim = c.numEntities(0);
  uint *idx  = new uint[3 * local_dim];
  uint *id  = new uint[3 * local_dim];
  real *tau_1_block = new real[3 * local_dim];  
  real *tau_2_block = new real[3 * local_dim];  
  real *normal_block = new real[3 * local_dim];

  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    ufc.update(*cell, mesh.distdata());
    
    (form.dofMaps())[1].tabulate_dofs(idx, ufc.cell, cell->index());
    
    uint ii = 0;
    uint jj = 0;    
    for(uint i = 0; i < 3; i++) 
    {
      for(VertexIterator v(*cell); !v.end(); ++v, ii++) 
      {
        if (!mesh.distdata().is_ghost(v->index(), 0)) 
        {
          tau_1_block[jj] = node_normal.tau_1[i].get(*v);
          tau_2_block[jj] = node_normal.tau_2[i].get(*v);
          normal_block[jj] = node_normal.normal[i].get(*v);
          id[jj++] = idx[ii];
        }
      }
    }

    tau_1.set(tau_1_block, jj, id);
    tau_2.set(tau_2_block, jj, id);
    normal.set(normal_block, jj, id);
  }

  tau_1.apply();
  tau_2.apply();
  normal.apply();
  delete[] tau_1_block;
  delete[] tau_2_block;
  delete[] normal_block;
  delete[] idx;
  delete[] id;
}

void NSESolver::merge(real *a,real *b,real *res,int an,int bn)
{
  real *ap,*bp,*rp;
  ap=a;
  bp=b;
  rp=res;

  while(ap<a+an && bp<b+bn){ 
    if(*ap <= *bp){
      *rp=*ap;
      ap++;
      rp++;
    }
    else { 
      *rp=*bp;
      rp++;
      bp++;
    }
  }
  if(ap<a+an){
    do
      *rp=*ap;
    while(++rp && ++ap<a+an);
  }
  else{
    do
      *rp=*bp;
    while(++rp && ++bp<b+bn);
  }
}

void NSESolver::ComputeLargestIndicators_cell(Mesh& mesh, Vector& e_indx,
        std::vector<int>& cells, real percentage)
{
  int N = mesh.numCells();
  int M = std::min((int)(N), 
                   (int)((real) 
                         (dolfin::MPI::numProcesses() > 1 ? 
                          mesh.distdata().global_numCells() : mesh.numCells()) * percentage * 0.01));
  
  if(dolfin::MPI::processNumber() == 0)
    dolfin_set("output destination","terminal");
  message("Computing largest indicators");
  message("percentage: %f", percentage);
  message("N: %d", N);
  message("M: %d", M);
  dolfin_set("output destination","silent");


  std::vector<std::pair<int, real> > indicators(N);
  real eind;
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    int id = (*cell).index();
    std::pair<int, real> p;
    p.first = id;
    uint ci = id;    
    if(dolfin::MPI::numProcesses() > 1)
      ci = mesh.distdata().get_cell_global(ci);
    e_indx.get(&eind, 1, &ci);      
    p.second = eind;    
    indicators[id] = p;
  }

  less_pair comp;
  std::sort(indicators.begin(), indicators.end(), comp);


  real *local_eind = new real[M];
  for(int i = 0; i < M; i++)
  {
    std::pair<int, real> p = indicators[N - 1 - i];
    local_eind[M - 1 - i] = p.second;
  }


  /*
   *  FIXME reduce memory usage
   *  merge only half of the received data
   */

  uint M_max, M_tot;
  MPI_Allreduce(&M, &M_max, 1, MPI_UNSIGNED, MPI_MAX, dolfin::MPI::DOLFIN_COMM);
  MPI_Allreduce(&M, &M_tot, 1, MPI_UNSIGNED, MPI_SUM, dolfin::MPI::DOLFIN_COMM);

  double *recv_eind = new double[M_max];
  double *global_eind = new double[M_tot];
  double *work = new double[M_tot];

  //  std::vector<double> global_eind;

  MPI_Status status;
  uint src,dest;
  uint rank =  dolfin::MPI::processNumber();
  uint size =  dolfin::MPI::numProcesses();
  uint nm = M;
  int num_recv;
  //  global_eind.insert(global_eind.begin(), local_eind, local_eind + M);
  std::memcpy(global_eind, local_eind, M*sizeof(real));

  for(uint i = 1; i < size; i++) {
    src =(rank - i + size) % size;
    dest = (rank + i) % size;

    MPI_Sendrecv(local_eind, M, MPI_DOUBLE, dest, 0, 
                 recv_eind, M_max, MPI_DOUBLE, src, 0, dolfin::MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_DOUBLE,&num_recv);
    //global_eind.insert(global_eind.end(), recv_eind, recv_eind + num_recv);
    merge(recv_eind, global_eind, work, num_recv, nm);
    std::memcpy(global_eind, work, M_tot * sizeof(real));
    nm += num_recv;
    
  }

  //  std::sort(global_eind.begin(), global_eind.end());
  cells.clear();
  int MM = (int)((real) (dolfin::MPI::numProcesses() > 1 ? 
                         mesh.distdata().global_numCells() : mesh.numCells()) * percentage * 0.01);
  int i = 0;
  for(int j = 0; j < MM; j++) {
    if( local_eind[M - 1 - i] >= global_eind[M_tot - 1 - j] ) {
      std::pair<int, real> p = indicators[N - 1 - i];
      cells.push_back(p.first);
      if( (i++) >= std::min(N, MM)) break;    
    }
  }

  dolfin_set("output destination", "terminal");
  message("%d marked cells on cpu %d", cells.size(), dolfin::MPI::processNumber());
  dolfin_set("output destination", "silent");

  
  delete[] local_eind;
  delete[] recv_eind;
  delete[] global_eind;
  delete[] work;
}

void NSESolver::ComputeRefinementMarkers(Mesh& mesh, real percentage,
        Vector& e_indx, MeshFunction<bool>& cell_refinement_marker)
{

  real error = 0.0;
  //ComputeError(error);

  //message("err: %g", error);
  
  std::vector<int> cells;
  ComputeLargestIndicators_cell(mesh, e_indx, cells, percentage);
    
  cell_refinement_marker.init(mesh, mesh.topology().dim());
  cell_refinement_marker = false;
    
  int M = cells.size();
       
  for(int i = 0; i < M; i++)
  {
    cell_refinement_marker.set(cells[i], true);
  }
}

// Utility function for mesh rotation
void NSESolver::computeX(Function& XX, Form* aM, Mesh& mesh)
{
  // Copy mesh coordinates into X array/function
  int d = mesh.topology().dim();
  UFC ufc(aM->form(), mesh, aM->dofMaps());
  Cell c(mesh, 0);
  uint local_dim = c.numEntities(0);
  uint *idx  = new uint[d * local_dim];
  uint *id  = new uint[d * local_dim];
  real *XX_block = new real[d * local_dim];

  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    ufc.update(*cell, mesh.distdata());
    (aM->dofMaps())[0].tabulate_dofs(idx, ufc.cell, cell->index());

    uint ii = 0;
    uint jj = 0;
    for(uint i = 0; i < d; i++)
    {
      for(VertexIterator v(*cell); !v.end(); ++v, ii++)
      {
        if (!mesh.distdata().is_ghost(v->index(), 0))
        {
          XX_block[jj] = v->x()[i];
          id[jj++] = idx[ii];
        }
      }
    }
    XX.vector().set(XX_block, jj, id);
  }
  XX.vector().apply();
  XX.sync_ghosts();
  delete[] XX_block;
  delete[] idx;
  delete[] id;
}

// Utility function for mesh rotation
void NSESolver::rotate(Function& XX, Form* aM, Mesh& mesh, double theta)
{
  MeshGeometry& geometry = mesh.geometry();

  uint d = mesh.topology().dim();
  uint N = mesh.numVertices();
  if(dolfin::MPI::numProcesses() > 1)
    N = mesh.distdata().global_numVertices();
  UFC ufc(aM->form(), mesh, aM->dofMaps());
  Cell c(mesh, 0);
  uint local_dim = c.numEntities(0);
  uint *idx  = new uint[d * local_dim];
  uint *id  = new uint[d * local_dim];
  real *XX_block = new real[d * local_dim];

  // Update the mesh                                                                                                                                                                                        
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    ufc.update(*cell, mesh.distdata());
    (aM->dofMaps())[0].tabulate_dofs(idx, ufc.cell, cell->index());

    XX.vector().get(XX_block, d * local_dim, idx);

    std::vector<double> xx, yy, zz;
    uint jj = 0;
    for(VertexIterator v(*cell); !v.end(); ++v)
    {
      for(unsigned int i = 0; i < d; i++)
      {
        if (i==0)
          xx.push_back(XX_block[i * local_dim + jj]);
        if (i==1)
          yy.push_back(XX_block[i * local_dim + jj]);
      }
      jj++;
    }
    uint j = 0;
    for(VertexIterator v(*cell); !v.end(); ++v)
    {
      Vertex& vertex = *v;

      real theta2;

      Point cp(0.35, -0.06, vertex.point()[2]);
      Point pdiff = vertex.point() - cp;

      real r = pdiff.norm();
      if(r < 0.5)
        theta2 = (0.0 - theta)*2*DOLFIN_PI/360.0;
      else if(r >= 0.5 && r < 1.0)
        theta2 = (0.0 - theta)*2*DOLFIN_PI/360.0*(1.0 - r) / (1.0 - 0.5);
      else
        theta2 = 0.0*2*DOLFIN_PI/360.0;

      for(unsigned int i = 0; i < d; i++)
      {
        if (i==0)
          XX_block[i * local_dim + j] = xx[j]*cos(theta2) - yy[j]*sin(theta2);
        if (i==1)
          XX_block[i * local_dim + j] = xx[j]*sin(theta2) + yy[j]*cos(theta2);

        geometry.x(vertex.index(), i) = XX_block[i * local_dim + j];
      }
      j++;
    }
  }

  delete[] XX_block;
  delete[] idx;
  delete[] id;

  MPI_Barrier(dolfin::MPI::DOLFIN_COMM);
}

NSESolver::~NSESolver ()
{
    delete assembler;
    delete h;
    delete n;
    delete nn;
    delete cv;
    
    bcs_m.clear();
    bcs_dm.clear();
    delete dbc_c;
    
    delete u;
    delete u0;
    delete p;
    delete p0;
    delete nuf;
    delete kf;
    delete c1f;
    delete c2f;
    delete hminf;
    delete dtu;
    
    delete rd_u;
    delete rd_p;
    
    delete up;
    delete pp;
    delete up0;
    delete dtup;
    
    delete Rm;
    delete Rmtot;
    delete Rc;
    delete Rctot;
    delete wm;
    delete wmtot;
    delete wc;
    delete wctot;
    
    delete ei_m;
    delete ei_c;
    delete eij_m;
    delete eij_c;
    delete eif;
    
    delete tau_1;
    delete tau_2;
    delete normal;
    
    delete ei;
    
    delete ap_m,
    delete Lp_m,
    
    delete ap_c;
    delete Lp_c,
    
    delete ad_m,
    delete Ld_m,
    
    delete ad_c,
    delete Ld_c,
    
    delete Lrep_m;
    delete Lrep_c;
    
    delete Md;
    delete Ml;
    
    delete MH1;
    delete MH12;
    delete MHerrest;
    delete MHerrestg;
    
    delete LRm;
    delete MRm;
    delete LRc;
    delete MRc;
    delete MRgm;
    delete MRgc;
    
    delete Lwm;
    delete Mwm;
    delete Lwc;
    delete Mwc;
    
    delete Mgwm;
    delete Mgwc;
    
    delete pdep_m;
    delete pdep_c;
    
    delete pded_m;
    delete pded_c;
    
    delete U;
    delete P;
    
    delete file_u;
    delete file_p;
    delete file_du;
    delete file_dp;
    delete file_m;
    
    delete Up;
    delete dtUp;
    delete Pp;
    delete Rmp;
    delete Rcp;
}


#warning "unused function"
void ComputeLargestIndicators_eind(Mesh& mesh, Vector& e_indx, std::vector<int>& cells,
                                                  real percentage)
{
  int N = mesh.numCells();
  real eind, sum_e, sum_e_local, max_e, max_e_local, min_e, min_e_local;
  sum_e = sum_e_local = max_e_local = 0.0;
  min_e_local = 1e6;
  
  std::vector<std::pair<int, real> > indicators(N);

  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    int id = (*cell).index();
    std::pair<int, real> p;
    p.first = id;
    uint ci = id;    
    if(dolfin::MPI::numProcesses() > 1)
      ci = mesh.distdata().get_cell_global(ci);
    e_indx.get(&eind, 1, &ci);      
    // Take absolute value
    eind = abs(eind);
    p.second = eind;    
    indicators[id] = p;
    max_e_local = std::max(max_e_local, eind);
    min_e_local = std::min(min_e_local, eind);
    sum_e_local += p.second;
  }

  less_pair comp;
  std::sort(indicators.begin(), indicators.end(), comp);

  MPI_Allreduce(&sum_e_local, &sum_e, 1, MPI_DOUBLE,
                MPI_SUM, dolfin::MPI::DOLFIN_COMM);

  MPI_Allreduce(&max_e_local, &max_e, 1, MPI_DOUBLE, 
                MPI_MAX, dolfin::MPI::DOLFIN_COMM);

  MPI_Allreduce(&min_e_local, &min_e, 1, MPI_DOUBLE, 
                MPI_MIN, dolfin::MPI::DOLFIN_COMM);

  real threshold = (percentage * 0.01 * sum_e);
  real cutoff = (max_e + min_e) / 2.0;
  real acc_local, acc;
  acc_local = acc = 0.0;

  int iter = 0;
  while ( (fabs(acc - threshold) / threshold )  > 1e-2  && (iter++) < 10)
  {
    cutoff = (max_e + min_e) / 2.0;
    acc = acc_local = 0.0;
    cells.clear();

    for (int i = 0; i < N; i++) 
    {
      std::pair<int, real> p = indicators[N - 1 - i];

      cells.push_back(p.first);
      acc_local += p.second;

      if ( p.second < cutoff )
        break;     
    }

    MPI_Allreduce(&acc_local, &acc, 1, MPI_DOUBLE, 
                  MPI_SUM, dolfin::MPI::DOLFIN_COMM);
        
    ( acc > threshold ? (min_e = cutoff ) : (max_e = cutoff));    
  }
}




#endif /* end of include guard: NSESOLVER_H */
