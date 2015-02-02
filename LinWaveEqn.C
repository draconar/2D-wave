//
// File:        Linear Wave Eqn model
// Description: SAMRAI interface class for linear wave eqn example.
//

#include "LinWaveEqn.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#include "tbox/Array.h"
#include "BoundaryBox.h"
#include "CartesianPatchGeometry.h"
#include "Index.h"
#include "tbox/IEEE.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "NodeData.h"
#include "NodeIndex.h"
#include "NodeIterator.h"
#include "Patch.h"
#include "PatchData.h"
#include "RefineOperator.h"
#include "tbox/Utilities.h"
#include "VariableDatabase.h"
#include "PatchNodeDataNormOpsReal.h"
//#define USE_FORTRAN
//#define USE_C

//CELLDATA
#include "HierarchyCellDataOpsReal.h"
#include "PatchCellDataOpsReal.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellVariable.h"

#define USE_CPP
//
using namespace SAMRAI;
using namespace tbox;
using namespace mesh;
using namespace geom;
using namespace appu;
using namespace hier;
using namespace xfer;
using namespace pdat;


//#define DEBUG_PRINT_ADVANCE
//#define DEBUG_PRINT_BOUNDARY
//#define DEBUG_PRINT_INITIAL

/*
 * Fortran function(s) used, if fortran is desired.
 */
#ifdef USE_FORTRAN
extern "C" {
#if (NDIM == 2)
void updateu_(const int&, const int&,
              const int&, const int&,
              const int&, const int&,
              const double*, const double&,
              double*,
              const double*, 
              const double*,
              const double*);
}
#endif
#endif


/*************************************************************************
 *
 * Constructor and Destructor.
 *
 ************************************************************************/

LinWaveEqn::LinWaveEqn(
   const string& object_name,
   Pointer<Database> input_db,
   Pointer<CartesianGridGeometry<NDIM> > grid_geom)
{

   /*
    * Set private data maintained by this class.  Note:  We use "d_" to 
    * denote PRIVATE data maintained by the class.  Variables that do
    * not contain "d_" are temporary variables used only in this method.
    */
   d_object_name = object_name;

   d_grid_geometry = grid_geom;   
	
  d_adaption_threshold = 0.2;
  d_scalar = 0.0;   

   /*
    * Initialize object with data read from input databases.
    */
   getFromInput(input_db);


}

LinWaveEqn::~LinWaveEqn()
{   
}


/*************************************************************************
 *
 * Setup variables used in the problem.
 *
 ************************************************************************/
void
LinWaveEqn::setupVariables(
   Pointer<VisItDataWriter<NDIM> > viz_writer)
{
   /* 
    * Set up variables and contexts.  The solution approach being used
    * (in 1D) is:
    *
    *     u^(n+1)_i = dt^2/h^2*(u^n_i+1 - 2u^n_i + u^n_i-1) 
    *                 + 2*u^n_i - u^(n-1)_i + dt^2*f
    *
    * so we need to store u at multiple time levels - n+1, n, n-1.
    * In SAMRAI, we refer to these as "contexts".  That is, we create
    * one node-centered variable "u", and three contexts:
    *    
    * NEW     - u^(n+1) (zero ghosts)
    * CURRENT - u^n     (1 ghost)
    * OLD     - u^(n-1) (zero ghosts)
    *
    * After each step of the above algorithm, we exchange NEW->CURRENT, 
    * CURRENT->OLD, OLD is thrown away.
    *
    * Below we construct u variable and its contexts.
    */
   VariableDatabase<NDIM>* variable_db = VariableDatabase<NDIM>::getDatabase();

   d_u_var = new NodeVariable<NDIM,double>("u",1);
   d_exact_var = new NodeVariable<NDIM,double>("exact",1);
   d_error_var = new NodeVariable<NDIM,double>("error",1);
   
   d_new_cxt = variable_db->getContext("NEW");
   d_cur_cxt = variable_db->getContext("CURRENT");
   d_old_cxt = variable_db->getContext("OLD");

   /*
    * We register variable-context pairs with SAMRAI's VariableDatabase<NDIM>.
    * Registering it identifies certain information about the variable,
    * namely its type (cell, node, face, etc.) and number of ghost cells
    * which will later be used by SAMRAI for communication operations,
    * so this is an important step. When a variable-context pair is 
    * registered, it returns an integer identifier that we can use later
    * to access data from the patch.
    */ 

   IntVector<NDIM> zero_ghosts(0);
   IntVector<NDIM> one_ghost(1);
   d_u_new_id = variable_db->registerVariableAndContext(d_u_var,
                                                        d_new_cxt,
                                                        zero_ghosts);
   d_u_cur_id = variable_db->registerVariableAndContext(d_u_var,
                                                        d_cur_cxt,
                                                        one_ghost);
   d_u_old_id = variable_db->registerVariableAndContext(d_u_var,
                                                        d_old_cxt,
                                                        zero_ghosts);
   d_exact_cur_id = variable_db->registerVariableAndContext(d_exact_var,
                                                        d_cur_cxt,
                                                        zero_ghosts);
   d_error_cur_id = variable_db->registerVariableAndContext(d_error_var,
                                                        d_cur_cxt,
                                                        zero_ghosts);

//tagging
//tbox::Pointer<hier::VariableContext> context_persistent_ptr(&d_context_persistent,false);
/*
d_scalar_persistent 
      = variable_db->registerVariableAndContext(
	Pointer< Variable<NDIM> >(&d_scalar,false)
	, context_persistent_ptr
	, IntVector<NDIM>(1) /* ghost cell width is 1 for stencil widths 
	);
*/
    /*
    * Register u values to be written by the viz writer.
    */
   viz_writer->registerPlotQuantity("U", "SCALAR",    d_u_cur_id );
   viz_writer->registerPlotQuantity("Exact", "SCALAR", d_exact_cur_id);
   viz_writer->registerPlotQuantity("Error", "SCALAR", d_error_cur_id);
}


/*************************************************************************
 *
 * Setup communication.
 *
 ************************************************************************/
void
LinWaveEqn::setupCommunication(
   const Pointer<PatchHierarchy<NDIM> > hierarchy) 
{
   /*
    * SAMRAI uses the notion of a grid-independent communication algorithm
    * to specify filling of ghost cells on a level, refining data from 
    * coarser levels to fill ghosts on a finer level, and coarsening data
    * to coarser levels.  A particular refine or coarsen algorithm
    * may be used throughout the calculation, while the grid is changing.
    *
    * Once the state of the grid hierarchy is known, the refine or coarsen
    * algorithm is used to construct a refine or coarsen schedule.  The 
    * schedule is effectively a specification of which patches talk to whom,
    * and may be used over and over again as long as the grid is static.
    * Once the grid configuration changes, the schedule must be reconstructed.
    *
    * Finally, the refine/coarsen schedule invokes the actual ghost fill
    * operation.  This is done inside the "advance()" method.
    */

   /*
    * Below, we construct a RefineAlgorithm<NDIM> to fill ghosts of neighboring
    * patches on a level and, ghosts of patches that live at a coarse
    * fine boundary. Note that the hierarchy is NOT needed for this operation,
    * so this step could be completed before the hierarchy is built.
    */ 

   Pointer<RefineOperator<NDIM> > refine_op = d_grid_geometry->
      lookupRefineOperator(d_u_var, "LINEAR_REFINE");

   Pointer<RefineAlgorithm<NDIM> > bdry_fill_alg = new RefineAlgorithm<NDIM>();
   bdry_fill_alg->registerRefine(d_u_cur_id,  // dest
                                 d_u_cur_id,  // source
                                 d_u_cur_id,  // scratch
                                 refine_op);

   /*
    * Now, build refine schedules for each level in the hierarchy.  Note
    * that the hierarchy IS needed for this operation.
    */
   d_bdry_fill_sched.resizeArray(hierarchy->getNumberLevels());
   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; ln--) {
      Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

      // Note:  a pointer to "this" tells the refine schedule to invoke
      // the setPhysicalBCs() method defined in this class.
      d_bdry_fill_sched[ln] = bdry_fill_alg->createSchedule(level,
                                                            ln-1,
                                                            hierarchy,
                                                            this);
   }
   
      
}

/*************************************************************************
 *
 * Advance hierarchy one timestep.
 *
 ************************************************************************/

void LinWaveEqn::advance(
   const Pointer<PatchHierarchy<NDIM> > hierarchy,
   const double dt,
   const double time)
{

   /*
    * First exchange ghost data.  Then compute u^(n+1)...
    *
    *     u^(n+1)_i = dt^2/h^2*(u^n_i+1 - 2u^n_i + u^n_i-1) 
    *                 + 2*u^n_i - u^(n-1)_i
    */

   // loop over hierarchy levels   
   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; ln--) {
      Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

      // fill ghosts of u_cur
      d_bdry_fill_sched[ln]->fillData(time);

      // allocate "new" data
      level->allocatePatchData(d_u_new_id, time);

      // loop over patches on level
      for (PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
         Pointer<Patch<NDIM> > patch = level->getPatch(ip());     

         // values for u - u_n+1, u_n, u_n-1         
         Pointer< NodeData<NDIM,double> > u_new = patch->getPatchData(d_u_new_id);
         Pointer< NodeData<NDIM,double> > u_cur = patch->getPatchData(d_u_cur_id);
         Pointer< NodeData<NDIM,double> > u_old = patch->getPatchData(d_u_old_id);

         // set values for f
         Box<NDIM> pbox = patch->getBox();
         const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = 
            patch->getPatchGeometry();
         const double* dx = patch_geom->getDx();
         const double* xdomainlo = d_grid_geometry->getXLower();
         double xloc[NDIM];
        
         IntVector<NDIM> zero_ghosts(0);
         NodeData<NDIM,double> f(pbox, 1, zero_ghosts); //verificar com cuidado
         for (NodeIterator<NDIM> i(pbox); i; i++) {
            NodeIndex<NDIM> node = i();  // i,j
            xloc[0] = xdomainlo[0] + dx[0]*((double)(node(0))); // x
            xloc[1] = xdomainlo[1] + dx[1]*((double)(node(1))); // y
            double pi = 3.141592654;
            double w = 2. * pi;
            (f)(node) = 0*w*w * sin(w*(xloc[0]-time)) * sin(w*xloc[1]);
         }

         // update u
         IntVector<NDIM> ghost_cells = u_cur->getGhostCellWidth();

#ifdef USE_FORTRAN
         const Index<NDIM> ifirst =patch->getBox().lower();
         const Index<NDIM> ilast  =patch->getBox().upper();

         updateu_(ifirst(0),ilast(0),
                  ifirst(1),ilast(1),
                  ghost_cells(0), ghost_cells(1),
                  dx,dt,
                  u_new->getPointer(),
                  u_cur->getPointer(),
                  u_old->getPointer(),
                  f.getPointer());         
#endif

#ifdef USE_C
         double* unew = u_new->getPointer();
         double* ucur = u_cur->getPointer();
         double* uold = u_old->getPointer();
         double* fdat = f.getPointer();

         const Index<NDIM> ifirst =pbox.lower();
         const Index<NDIM> ilast  =pbox.upper();

         Index<NDIM> nxg = ilast - ifirst + IntVector<NDIM>(2)  // # nodes in X,Y
            + ghost_cells + ghost_cells;            //   + ghost cells

         Index<NDIM> nx = ilast - ifirst + IntVector<NDIM>(2);  // # nodes in X,Y
         /*
         for (int ic1 = 0; ic1 < nx[1]; ic1++) {
            for (int ic0 = 0; ic0 < nx[0]; ic0++) {
               int ij_nogh = ic1*nx[0] + ic0;
               int ij_withg = (ic1+ghost_cells(1))*nxg[0] 
                  + (ic0+ghost_cells(0));
               double xdiff  = (ucur[ij_withg+1] - 2.0*ucur[ij_withg]
                  + ucur[ij_withg-1])/(dx[0]*dx[0]);
               double ydiff  = (ucur[ij_withg+nxg[0]] - 2.0*ucur[ij_withg]
                  + ucur[ij_withg-nxg[0]])/(dx[1]*dx[1]);
               
               unew[ij_nogh] = dt*dt*(xdiff + ydiff) 
                  + 2.0*ucur[ij_withg] - uold[ij_nogh] + dt*dt*fdat[ij_nogh];
               
            }
         }*/
#endif
         
#ifdef USE_CPP
//METHOD SOLVER IS HERE
         for (NodeIterator<NDIM> i(pbox); i; i++) {
            NodeIndex<NDIM> node = i();  // i,j
            NodeIndex<NDIM> ip1_x = node;
            NodeIndex<NDIM> im1_x = node;
            NodeIndex<NDIM> ip1_y = node;
            NodeIndex<NDIM> im1_y = node;
            ip1_x(0) += 1;  // i+1,j
            im1_x(0) -= 1;  // i-1,j
            ip1_y(1) += 1;  // i,j+1
            im1_y(1) -= 1;  // i,j-1

            double xdiff = (dt/dx[0]) * ( (*u_cur)(im1_x) - (*u_cur)(ip1_x) ) 
                            + (dt/dx[0]) * (dt/dx[0]) * ( (*u_cur)(im1_x) - 2 * (*u_cur)(node) + (*u_cur)(ip1_x) );

            double ydiff = (dt/dx[1]) * ( (*u_cur)(im1_y) - (*u_cur)(ip1_y) ) 
                            + (dt/dx[1]) * (dt/dx[1]) * ( (*u_cur)(im1_y) - 2 * (*u_cur)(node) + (*u_cur)(ip1_y) );

	    (*u_new)(node) = (*u_cur)(node) + (xdiff + ydiff)/2;
            //(*u_new)(node) = (xdiff + ydiff - (*u_old)(node) - dt * (*u_old)(node) - 0.5 * dt * dt * (*u_old)(node) + 0.5 *dt*dt* (*u_cur)(node) ) / ( dt * (0.5 * dt - 1) );
	}
#endif            
         
      } // end loop over patches

      /*
       * Update state of variables by switching contexts:
       *   1) transfer "cur" to "old"
       *   2) transfer "new" to "cur"
       */
      for (PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {

         Pointer<Patch<NDIM> > patch = level->getPatch(ip());  
         Pointer< NodeData<NDIM,double> > u_new = patch->getPatchData(d_u_new_id);
         Pointer< NodeData<NDIM,double> > u_cur = patch->getPatchData(d_u_cur_id);
         Pointer< NodeData<NDIM,double> > u_old = patch->getPatchData(d_u_old_id);

         u_old->copy(*u_cur);
         u_cur->copy(*u_new);

#ifdef DEBUG_PRINT_ADVANCE
         plog << "----Advance - cur"  << endl;
         plog << "Patch<NDIM>: " << patch->getPatchNumber() 
              << "  Box<NDIM>: " << patch->getBox() << endl;
         u_cur->print(u_cur->getGhostBox());
         
         plog << "----Advance - old"  << endl;
         plog << "Patch<NDIM>: " << patch->getPatchNumber() 
              << "  Box<NDIM>: " << patch->getBox() << endl;
         u_old->print(u_old->getGhostBox());
#endif
 
      }

      /*
       * We no longer need the "new" context so we de-allocate it.
       */
      level->deallocatePatchData(d_u_new_id);

   } // loop over levels
   
}


/*************************************************************************
 *
 * Determine timestep to be used - CFL = k/h => k = CFL*h
 *
 ************************************************************************/
double
LinWaveEqn::getDt( const Pointer<PatchHierarchy<NDIM> > hierarchy)
{
   double dx_min = 1.e12;
   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; ln--) {
      Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
      for (PatchLevel<NDIM>::Iterator p(level); p; p++) {
         Pointer<Patch<NDIM> > patch = level->getPatch(p());
         
         Box<NDIM> pbox = patch->getBox();
         const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = 
            patch->getPatchGeometry();
         const double* dx = patch_geom->getDx();
         
         for (int i = 0; i < NDIM; i++) {
	   if(dx[i] < dx_min) dx_min = dx[i];//Utilities::dmin(dx_min,dx[i]);
         }
      }
   }

   return(dx_min*d_cfl);
}

/*************************************************************************
 *
 * Compare solution to exact solution
 *
 ************************************************************************/
void
LinWaveEqn::compareToExact(  const Pointer<PatchHierarchy<NDIM> > hierarchy,
                             const double time)
{
   /*
    * Exact solution is:
    *
    *   u(x,y,t) = sin((x-t) * (y-t))
    *
    * Compute the RMS of this solution...
    */
  for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; ln--) {
     int ln =0;
     Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

      // fill ghosts of u_cur
      d_bdry_fill_sched[ln]->fillData(time);

     int nodes_on_level = 0;
     double rms_on_level = 0.;
     double max_norm_on_level = 0.;
     NodeIndex<NDIM> max_loc;
     for (PatchLevel<NDIM>::Iterator p(level); p; p++) {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());

        Box<NDIM> pbox = patch->getBox();
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = 
           patch->getPatchGeometry();
        const double* dx = patch_geom->getDx();
        const double* xdomainlo = d_grid_geometry->getXLower();
        const Index<NDIM> ifirst = pbox.lower();
        const Index<NDIM> ilast  = pbox.upper();
        double xloc[NDIM];
        
        /*
         * Access u_n, exact, error data from patch.
         */
        Pointer< NodeData<NDIM,double> > u_cur = 
           patch->getPatchData(d_u_cur_id);
        Pointer< NodeData<NDIM,double> > exact = 
           patch->getPatchData(d_exact_cur_id);
        Pointer< NodeData<NDIM,double> > error = 
           patch->getPatchData(d_error_cur_id);

        Box<NDIM> ghost_box = u_cur->getGhostBox();

        double pi = 3.141592654;
        double w = 4. * pi;
        double exact_soln = 0.;
        double diff = 0.;


        for (NodeIterator<NDIM> i(pbox); i; i++) {
           NodeIndex<NDIM> node = i();  // i,j
           xloc[0] = xdomainlo[0] + dx[0]*((double)(node(0))); // x
           xloc[1] = xdomainlo[1] + dx[1]*((double)(node(1))); // y
           //exact_soln = sin(w*(xloc[0]-time)); // 1D
           exact_soln = sin( w * (xloc[0]-time) ); // 2D
           diff       = (*u_cur)(node) - exact_soln;
           (*exact)(node) = exact_soln;
           (*error)(node) = diff;
           rms_on_level = rms_on_level + diff*diff;
           nodes_on_level += 1;
           if (diff > max_norm_on_level) {
              max_norm_on_level = diff;
              max_loc = node;
           }
        }

     }

     //pout.precision(12);

//BEGIN SETUP L2
/*math::PatchNodeDataNormOpsReal<NDIM, double> *l2 = new math::PatchNodeDataNormOpsReal<NDIM,double>;
double l2result;

	l2result = l2->L2Norm(u_cur, pbox, NULL);
	pout << "L2 == " << l2result << endl;
*/
//END NORM2          

	pout << "RMS of error on level " << level->getLevelNumber()
          << " at time " << time 
          << ": " << rms_on_level/(double)nodes_on_level << endl;
     pout << "MAX norm of error: " << max_norm_on_level
          << " at location " << max_loc
          << endl;

  }

}

      

/*************************************************************************
 *
 * Methods inherited from StandardTagAndInitStrategy<NDIM>.
 *
 ************************************************************************/


/*
 * Allocate storage and initialize data on the level.
 */
void
LinWaveEqn::initializeLevelData(
   const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
   const int level_number,
   const double time,
   const bool can_be_refined,
   const bool initial_time,
   const Pointer<BasePatchLevel<NDIM> > old_level,
   const bool allocate_data)
{

   /*
    * Set initial data on hierarchy level.  We will allocate data for
    * all three contexts of u - new, curr, old - and initialize
    * u^n (u_cur) and u^n-1 (u_old).
    */
   Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
   Pointer<PatchLevel<NDIM> > level =
      hierarchy->getPatchLevel(level_number);

   /*
    * Allocate storage for "cur" and "old" data.  
    */
   if (allocate_data) {
      level->allocatePatchData(d_u_cur_id, time);
      level->allocatePatchData(d_u_old_id, time);
      level->allocatePatchData(d_exact_cur_id, time);
      level->allocatePatchData(d_error_cur_id, time);
   }
   
   /*
    * Set initial state on all patches of the level.
    *
    *   u(x,y,t) = sin(w(x-t))sin(wy) where w = 4*pi
    *
    * At t=0,
    *
    *   u(x,y,0) = sin(wx)sin(wy) = f(x,y)
    *   du/dt(x,y,0) = -w cos(wx)sin(wy) = g(x,y)
    */
   for (PatchLevel<NDIM>::Iterator p(level); p; p++) {
      Pointer<Patch<NDIM> > patch = level->getPatch(p());

      Box<NDIM> pbox = patch->getBox();
      const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = 
         patch->getPatchGeometry();
      const double* dx = patch_geom->getDx();
      const double* xdomainlo = d_grid_geometry->getXLower();
      const Index<NDIM> ifirst = pbox.lower();
      const Index<NDIM> ilast  = pbox.upper();
      double xloc[NDIM];

      /*
       * Set u_n data - u_n = f(x,y)
       */
      Pointer< NodeData<NDIM,double> > u_cur = 
         patch->getPatchData(d_u_cur_id);

      u_cur->fillAll(0.);
      
      double pi = 3.141592654;
      double w = 4. * pi;

      for (NodeIterator<NDIM> i(pbox); i; i++) {
         NodeIndex<NDIM> node = i();  // i,j
         xloc[0] = xdomainlo[0] + dx[0]*((double)(node(0))); // x
         xloc[1] = xdomainlo[1] + dx[1]*((double)(node(1))); // y
         //pout << "x[" << node << "]:  " << xloc[0] << endl;
         //(*u_cur)(node) = sin(w*xloc[0]); // 1D
	if ( (xloc[1] > 0.5 && xloc[1] < 0.8) && (xloc[0] > 0.5 && xloc[0] < 0.8)) {
         (*u_cur)(node) = 1+0*sin( w*xloc[0] ); // 2D
	} else {
	 (*u_cur)(node) = 0*sin( w*xloc[0] ); // 2D	
	}	
      }

      /*
       * Set u_(n-1) data 
       */
      Pointer< NodeData<NDIM,double> > u_old = 
         patch->getPatchData(d_u_old_id);
      
      double dudt = 0.;
      double dt = getDt(hierarchy);
      for (NodeIterator<NDIM> i(pbox); i; i++) {
         NodeIndex<NDIM> node = i();  // i,j
         xloc[0] = xdomainlo[0] + dx[0]*((double)(node(0))); // x
         xloc[1] = xdomainlo[1] + dx[1]*((double)(node(1))); // y
         //(*u_old)(node) = sin(w*(xloc[0]+dt)); // 1D
	if( (xloc[1] > 0.5 && xloc[1] < 0.8) && (xloc[0] > 0.5 && xloc[0] < 0.8)){
         (*u_cur)(node) = 1+0*sin( w*xloc[0] ); // 2D
	} else {
	 (*u_cur)(node) = 0*sin( w*xloc[0] ); // 2D	
	}
      }

#ifdef DEBUG_PRINT_INITIAL
      plog << "----Initialize level data - cur"  << endl;
      plog << "Patch<NDIM>: " << patch->getPatchNumber() 
           << "  Box<NDIM>: " << patch->getBox() << endl;
      u_cur->print(u_cur->getGhostBox());

      plog << "----Initialize level data - old"  << endl;
      plog << "Patch<NDIM>: " << patch->getPatchNumber() 
           << "  Box<NDIM>: " << patch->getBox() << endl;
      u_old->print(u_old->getGhostBox());
#endif
      
   }
   
}

/*
 * Perform operations necessary when grid changes for dynamic grid
 * calculations (to be added later...)
 */
void 
LinWaveEqn::resetHierarchyConfiguration(
   const Pointer<BasePatchHierarchy<NDIM> >   hierarchy,
   const int coarsest_level,
   const int finest_level)
{
   (void) hierarchy;
   (void) coarsest_level;
   (void) finest_level;
}

/*
 * Routine to do cell tagging.  Add later...
 */
void
LinWaveEqn::applyGradientDetector(
   const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
   const int level_number,
   const double time,
   const int tag_index,
   const bool initial_time,
   const bool uses_richardson_extrapolation_too)
{
   (void) hierarchy;
   (void) level_number;
   (void) time;
   (void) tag_index;
   (void) initial_time;
   (void) uses_richardson_extrapolation_too;
}

/*
*************************************************************************
*
* Methods inherited from RefinePatchStrategy<NDIM>.
*
***********************************************************************
*/

/*
 * Set boundary conditions.
 */
void 
LinWaveEqn::setPhysicalBoundaryConditions(
   Patch<NDIM>& patch,
   const double time,
   const IntVector<NDIM>& ghost_width_to_fill)
{

   (void) time;

   const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
//   const double* dx = patch_geom->getDx();
//   const double* xpatchhi = patch_geom->getXUpper();
//   const double* xdomainhi = d_grid_geometry->getXUpper();

   /*
    * Get node and edge boundary boxes.
    */   
   const Array<BoundaryBox<NDIM> > node_bdry =
      patch_geom->getCodimensionBoundary(NDIM);
   const int num_node_bdry_boxes = node_bdry.getSize();

   const Array<BoundaryBox<NDIM> > edge_bdry =
      patch_geom->getCodimensionBoundary(NDIM-1);
   const int num_edge_bdry_boxes = edge_bdry.getSize();

   /*
    * Pointer to data in ghost regions.
    */
   Pointer< NodeData<NDIM,double> > u_cur = patch.getPatchData(d_u_cur_id);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!u_cur.isNull());
#endif
   IntVector<NDIM> ghost_cells = u_cur->getGhostCellWidth();
   const Box<NDIM> interior(patch.getBox());
   const double* dx = patch_geom->getDx();
   const double* xdomainlo = d_grid_geometry->getXLower();
   const Index<NDIM> ifirst = interior.lower();
   const Index<NDIM> ilast  = interior.upper();
   double xloc[NDIM];

   int i;
   double pi = 3.141592654;
   double w = 4. * pi;

   /*
    *  Set boundary conditions:
    *     u(x,y,t) = sin(xy)
    */   

   /*
    * Patch<NDIM> edges (i.e. 4 sides)
    */

   for ( i = 0; i < num_edge_bdry_boxes; i++ ) {
      /*
       * location index:
       *    0,1 - X lower,upper
       *    2,3 - Y lower,upper
       */    
      int loc_indx = edge_bdry[i].getLocationIndex();

      Box<NDIM> fill_box = patch_geom->getBoundaryFillBox(edge_bdry[i],
                                                    interior,
                                                    ghost_width_to_fill);

      /*
       * Iterate over edges of fill box.
       */
      for (NodeIterator<NDIM> i(fill_box); i; i++) {
         NodeIndex<NDIM> node = i();  //i,j

         switch(loc_indx) {
            case 0: {
               // X lower
               NodeIndex<NDIM> ll = NodeIndex<NDIM>(node, NodeIndex<NDIM>::LowerLeft); 
               xloc[0] = xdomainlo[0] + dx[0]*((double)(ll(0))); // x
               xloc[1] = xdomainlo[1] + dx[1]*((double)(ll(1))); // y
               //(*u_cur)(ll) = sin(w*(xloc[0]-time)); // 1D
               (*u_cur)(ll) = 0+0*sin( w*(xloc[0]-time)); // 2D
               NodeIndex<NDIM> ul = NodeIndex<NDIM>(node, NodeIndex<NDIM>::UpperLeft); 
               xloc[0] = xdomainlo[0] + dx[0]*((double)(ul(0))); // x
               xloc[1] = xdomainlo[1] + dx[1]*((double)(ul(1))); // y
               //(*u_cur)(ul) = sin(w*(xloc[0]-time)); // 1D         
               (*u_cur)(ul) = 0+0*sin( w*(xloc[0]-time)); // 2D
               //(*u_cur)(im1_x) = d_dirichlet_edge[1]; 
               break;
            }
            case 1: {
               // X upper
               NodeIndex<NDIM> lr = NodeIndex<NDIM>(node, NodeIndex<NDIM>::LowerRight); 
               xloc[0] = xdomainlo[0] + dx[0]*((double)(lr(0))); // x
               xloc[1] = xdomainlo[1] + dx[1]*((double)(lr(1))); // y
               //(*u_cur)(lr) = sin(w*(xloc[0]-time)); // 1D
               //(*u_cur)(lr) = sin( xloc[0] * xloc[1]); // 2D
               (*u_cur)(lr) = 0+0*sin( w*(xloc[0] -time)); // 2D
               NodeIndex<NDIM> ur = NodeIndex<NDIM>(node, NodeIndex<NDIM>::UpperRight); 
               xloc[0] = xdomainlo[0] + dx[0]*((double)(ur(0))); // x
               xloc[1] = xdomainlo[1] + dx[1]*((double)(ur(1))); // y
               //(*u_cur)(ur) = sin(w*(xloc[0]-time)); // 1D    
               //(*u_cur)(ur) = sin( xloc[0] * xloc[1]); // 2Dkta
               (*u_cur)(ur) = 0+0*sin(w*( xloc[0] - time)); // 2D
               //(*u_cur)(ip1_x) = d_dirichlet_edge[1];            
               break;
            }
            case 2: {
               // Y lower
               NodeIndex<NDIM> ll = NodeIndex<NDIM>(node, NodeIndex<NDIM>::LowerLeft); 
               xloc[0] = xdomainlo[0] + dx[0]*((double)(ll(0))); // x
               xloc[1] = xdomainlo[1] + dx[1]*((double)(ll(1))); // y
               //(*u_cur)(ll) = sin(w*(xloc[0]-time)); // 1D
               //(*u_cur)(ll) = sin( xloc[0] *xloc[1]); // 2D  Kta
               (*u_cur)(ll) = 0+0*sin(w*( xloc[0] -time)); // 2D
               NodeIndex<NDIM> lr = NodeIndex<NDIM>(node, NodeIndex<NDIM>::LowerRight); 
               xloc[0] = xdomainlo[0] + dx[0]*((double)(lr(0))); // x
               xloc[1] = xdomainlo[1] + dx[1]*((double)(lr(1))); // y
               //(*u_cur)(lr) = sin(w*(xloc[0]-time)); // 1D    
               //(*u_cur)(lr) = sin( xloc[0] * xloc[1]); // 2D
               (*u_cur)(lr) = 0+0*sin( w*(xloc[0]-time)); // 2D
               //(*u_cur)(im1_y) = d_dirichlet_edge[2];            
               break;
            }
            case 3: {
               // Y upper
               NodeIndex<NDIM> ul = NodeIndex<NDIM>(node, NodeIndex<NDIM>::UpperLeft);
               xloc[0] = xdomainlo[0] + dx[0]*((double)(ul(0))); // x
               xloc[1] = xdomainlo[1] + dx[1]*((double)(ul(1))); // y
               //(*u_cur)(ul) = sin(w*(xloc[0]-time)); // 1D
               (*u_cur)(ul) = 0+sin( w*xloc[0]); // 2D
               NodeIndex<NDIM> ur = NodeIndex<NDIM>(node, NodeIndex<NDIM>::UpperRight);
               xloc[0] = xdomainlo[0] + dx[0]*((double)(ur(0))); // x
               xloc[1] = xdomainlo[1] + dx[1]*((double)(ur(1))); // y
               //(*u_cur)(ur) = sin(w*(xloc[0]-time)); // 1D    
               (*u_cur)(ur) = 0+0*sin( w*(xloc[0] - time)); // 2D
               //(*u_cur)(ip1_y) = d_dirichlet_edge[3];            
               break;
            }
         }
      }
   } // loop over edge bdry boxes
      
#ifdef DEBUG_PRINT_BOUNDARY
    plog << "----Boundary Conditions "  << endl;
    plog << "Patch<NDIM>: " << patch.getPatchNumber() 
         << "  Box<NDIM>: " << patch.getBox() << endl;
    u_cur->print(u_cur->getGhostBox());
#endif

}




/*
*************************************************************************
*
* Get data from input database.                                         *
*                                                                       *
*************************************************************************
*/ 
void 
LinWaveEqn::getFromInput(Pointer<Database> input_db)
{

   d_cfl = input_db->getDoubleWithDefault("cfl",0.5);
   
//   if (d_cfl > 1/sqrt(2.)) {
//      TBOX_ERROR("CFL should be < 1/sqrt(2) - 0.707");
//   }
   
   d_initial_a = input_db->getDoubleWithDefault("initial_a",0.0);
   d_initial_b = input_db->getDoubleWithDefault("initial_b",1.0);

   if (input_db->keyExists("Boundary_data")) {
      Pointer<Database> boundary_db = input_db->getDatabase("Boundary_data");
      
      string bdry_type = "Dirichlet";
      bdry_type = boundary_db->getStringWithDefault("type", bdry_type);

      if (bdry_type == "Dirichlet") {

         d_dirichlet_edge[0] = boundary_db->
            getDoubleWithDefault("Xlo_edge_dirichlet",d_dirichlet_edge[0]);
         d_dirichlet_edge[1] = boundary_db->
            getDoubleWithDefault("Xhi_edge_dirichlet",d_dirichlet_edge[1]);
         d_dirichlet_edge[2] = boundary_db->
            getDoubleWithDefault("Ylo_edge_dirichlet",d_dirichlet_edge[2]);
         d_dirichlet_edge[3] = boundary_db->
            getDoubleWithDefault("Yhi_edge_dirichlet",d_dirichlet_edge[3]);

         d_dirichlet_node[0] = boundary_db->
            getDoubleWithDefault("XloXlo_corner_dirichlet",
                                 d_dirichlet_node[0]);
         d_dirichlet_node[1] = boundary_db->
            getDoubleWithDefault("XhiYlo_corner_dirichlet",
                                 d_dirichlet_node[1]);
         d_dirichlet_node[2] = boundary_db->
            getDoubleWithDefault("XloYhi_corner_dirichlet",
                                 d_dirichlet_node[2]);
         d_dirichlet_node[3] = boundary_db->
            getDoubleWithDefault("XhiYhi_corner_dirichlet",
                                 d_dirichlet_node[3]);

      } else {
         TBOX_ERROR("Boundary type " << bdry_type << " not defined.");
      }
   }   
      
}

/*
*************************************************************************
*                                                                       *
* Prints class data - writes out info in class if exception is thrown   *
*                                                                       *
*************************************************************************
*/

void LinWaveEqn::printClassData(ostream &os) const
{
   fflush(stdout);
   int j;

   os << "ptr LinWaveEqn = " << (LinWaveEqn*) this << endl;

   os << "d_object_name = " << d_object_name << endl;

   os << "d_cfl = " << d_cfl << endl;
   
   os << "d_u_new_id = " << d_u_new_id << endl;
   os << "d_u_cur_id = " << d_u_cur_id << endl;
   os << "d_u_old_id = " << d_u_old_id << endl;

   os << "d_initial_a = " << d_initial_a << endl;
   os << "d_initial_b = " << d_initial_b << endl;
   
   os << "Boundary Condition data..." << endl;
   for (j = 0; j < 2*NDIM; j++) {
      os << "j: " << j 
         << "  d_dirichlet_edge[j] = " << d_dirichlet_edge[j]
         << "\n  d_dirichlet_node[j] = " << d_dirichlet_node[j]
         << endl;
   }
   
}

