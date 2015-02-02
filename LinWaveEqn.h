//
// File:        LinWaveEqn.h
// Package:     SAMRAI mesh
// Copyright:   (c) 1997-2002 The Regents of the University of California
// Release:     $Name:  $
// Revision:    $Revision: 1.2 $
// Modified:    $Date: 2002/08/21 01:31:51 $
// Description: Model class for sample 2nd order linear wave eqn example.
//
 
#include "SAMRAI_config.h"

/*
 * Header file for base classes.
 */
#include "StandardTagAndInitStrategy.h"
#include "RefinePatchStrategy.h"
#include "CoarsenPatchStrategy.h"

/*
 * Header file for SAMRAI classes referenced in this class.
 */

#include "tbox/Array.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "tbox/Database.h"
#include "IntVector.h"
#include "tbox/IOStream.h"
#include "NodeVariable.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "VariableContext.h"
#include "RefineAlgorithm.h"
#include "RefineSchedule.h"
#include "VisItDataWriter.h"
#include <string.h>
#include <stdio.h>

//CELLDATA
#include "HierarchyCellDataOpsReal.h"
#include "PatchCellDataOpsReal.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellVariable.h"


using namespace std;
using namespace SAMRAI;
using namespace tbox;
using namespace mesh;
using namespace geom;
using namespace appu;
using namespace hier;
using namespace xfer;
using namespace pdat;


#ifndef NULL
#define NULL (0)
#endif

/*! \mainpage
 * The LinWaveEqn class provides an implementation that uses SAMRAI
 * classes to solve the 2nd order linear wave eqn (shown below in 1d):
 *
 *    \f$ u_{tt} - u_{xx} = f \f$ 
 * 
 * with initial conditions
 *
 *    \f$ u(x,0) = a(x) \f$
 *    \f$ u_{t}(x,0) = b(x) \f$
 * 
 * The class provides interfaces to problem dependent operations that
 * are expected by gridding operations performed in SAMRAI.  For example, 
 * initializing data on a level and resetting hierarchy after regrid and 
 * pre and post process operations performed when transferring data 
 * between hierarchy levels.  The class also provides an "advance()" method
 * which computes u^(n+1) on hierarchy levels.  
 *
 * Other methods may be added as needed.
 * 
 * Input Parameters:
 * 
 *    - \b cfl         
 *       double specifying CFL (k/h = wave freq/dx) for the problem.
 *
 *
 * A sample input entry might look like:
 *  
 *    LinWaveEqn {
 *
 *       cfl = 0.1
 *
 *    }
 */

/*!
 * classe principal
 */
class LinWaveEqn : 
   public StandardTagAndInitStrategy<NDIM>,
   public RefinePatchStrategy<NDIM>,
   public CoarsenPatchStrategy<NDIM>
{
public:
   /*!*
    * Default constructor for LinWaveEqn.
    */     
   LinWaveEqn(const string& object_name,
              Pointer<Database> input_db,
              Pointer<CartesianGridGeometry<NDIM> > grid_geom);

   /*!*
    * Empty destructor for LinWaveEqn.
    */
   virtual ~LinWaveEqn();

/*************************************************************************
 *
 * Methods particular to LinWaveEqn class.
 *
 ************************************************************************/

   /*!*
    * Setup variables used in the problem.
    */
   void setupVariables(
      Pointer<VisItDataWriter<NDIM> > viz_writer);   

   /*!*
    * Setup communication schedules for the functions on different levels.
    */
   void setupCommunication(const Pointer<PatchHierarchy<NDIM> > hierarchy);

   /*!*
    * Advance data one timestep on hierarchy.
    */
   void advance(const Pointer<PatchHierarchy<NDIM> > hierarchy,
                const double dt,
                const double time = 0.);

   /*!*
    * Determine timestep to be used for the grid hierarchy.
    */
   double getDt(const Pointer<PatchHierarchy<NDIM> > hierarchy);

   /*!*
    * Compare to exact solution.
    */
   void compareToExact(const Pointer<PatchHierarchy<NDIM> > hierarchy,
                       const double time);
   
   /*!*
    * Prints all class data members, if exception is thrown.
    */
   void printClassData(ostream &os) const;


/*!************************************************************************
 *
 * Methods inherited from StandardTagAndInitStrategy<NDIM>.
 *
 ************************************************************************/
 
   /**
    * Initialize data on a new level after it is inserted into an AMR patch
    * hierarchy by the gridding algorithm.  The level number indicates
    * that of the new level.
    *
    * Generally, when data is set, it is interpolated from coarser levels
    * in the hierarchy.  If the old level pointer in the argument list is
    * non-null, then data is copied from the old level to the new level
    * on regions of intersection between those levels before interpolation
    * occurs.   In this case, the level number must match that of the old 
    * level.  The specific operations that occur when initializing level 
    * data are determined by the particular solution methods in use; i.e.,
    * in the subclass of this abstract base class.
    *
    * The boolean argument initial_time indicates whether the level is
    * being introduced for the first time (i.e., at initialization time),
    * or after some regrid process during the calculation beyond the initial
    * hierarchy construction.  This information is provided since the
    * initialization of the data may be different in each of those
    * circumstances.  The can_be_refined boolean argument indicates whether
    * the level is the finest allowable level in the hierarchy.
    */
    
   virtual void
   initializeLevelData(const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                       const int level_number,
                       const double time,
                       const bool can_be_refined,
                       const bool initial_time,
                       const Pointer<BasePatchLevel<NDIM> > old_level = 
                             Pointer<BasePatchLevel<NDIM> >(NULL),
                       const bool allocate_data = true);

   /*!*
    * After hierarchy levels have changed and data has been initialized on 
    * the new levels, this routine can be used to reset any information 
    * needed by the solution method that is particular to the hierarchy 
    * configuration.  For example, the solution procedure may cache 
    * communication schedules to amortize the cost of data movement on the 
    * AMR patch hierarchy.  This function will be called by the gridding 
    * algorithm after the initialization occurs so that the algorithm-specific
    * subclass can reset such things.  Also, if the solution method must 
    * make the solution consistent across multiple levels after the hierarchy 
    * is changed, this process may be invoked by this routine.  Of course the 
    * details of these processes are determined by the particular solution 
    * methods in use.
    *
    * The level number arguments indicate the coarsest and finest levels
    * in the current hierarchy configuration that have changed.  It should
    * be assumed that all intermediate levels have changed as well.
    */
   virtual void resetHierarchyConfiguration(
      const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
      const int coarsest_level,
      const int finest_level);

   /*!*
    * Set tags to the specified tag value where refinement of the given
    * level should occur using the user-supplied gradient detector.  The 
    * value "tag_index" is the index of the cell-centered integer tag
    * array on each patch in the hierarchy.  The boolean argument indicates 
    * whether cells are being tagged on the level for the first time; 
    * i.e., when the hierarchy is initially constructed.  If it is false, 
    * it should be assumed that cells are being tagged at some later time 
    * after the patch hierarchy was initially constructed.  This information
    * is provided since the application of the error estimator may be 
    * different in each of those circumstances.
    */
   virtual void 
   applyGradientDetector(const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                         const int level_number,
                         const double time,
                         const int tag_index,
                         const bool initial_time,
                         const bool uses_richardson_extrapolation_too);

/*************************************************************************
 *
 * Methods inherited from RefinePatchStrategy<NDIM>.
 *
 ************************************************************************/
 
   /*!*
    * Set the data at patch boundaries corresponding to the physical domain
    * boundary.  The specific boundary conditions are determined by the user.
    */
   virtual void setPhysicalBoundaryConditions(
      Patch<NDIM>& patch,
      const double time,
      const IntVector<NDIM>& ghost_width_to_fill);

   /*!*
    * Perform user-defined refining operations.  This member function
    * is called before the other refining operators.  The preprocess
    * function should refine data from the scratch components of the
    * coarse patch into the scratch components of the fine patch on the
    * specified fine box region.  This version of the preprocess function
    * operates on a a single box at a time.  The user must define this
    * routine in the subclass.
    */
   virtual void preprocessRefine(
      Patch<NDIM>& fine,
      const Patch<NDIM>& coarse,
      const Box<NDIM>& fine_box,
      const IntVector<NDIM>& ratio)
      {
         (void) fine;
         (void) coarse;
         (void) fine_box;
         (void) ratio;
      }

   /*!*
    * Perform user-defined refining operations.  This member function
    * is called after the other refining operators.  The postprocess
    * function should refine data from the scratch components of the
    * coarse patch into the scratch components of the fine patch on the
    * specified fine box region.  This version of the postprocess function
    * operates on a a single box at a time.  The user must define this
    * routine in the subclass.
    */
   virtual void postprocessRefine(
      Patch<NDIM>& fine,
      const Patch<NDIM>& coarse,
      const Box<NDIM>& fine_box,
      const IntVector<NDIM>& ratio)
      {
         (void) fine;
         (void) coarse;
         (void) fine_box;
         (void) ratio;
      }

   /*!*
    * Return maximum stencil width needed for user-defined
    * data interpolation operations.  Default is to return
    * zero, assuming no user-defined operations provided.
    */
   virtual IntVector<NDIM> getRefineOpStencilWidth() const
   {
      return(IntVector<NDIM>(0));
   }

/*!************************************************************************
 *
 * Methods inherited from CoarsenPatchStrategy<NDIM>.
 *
 ************************************************************************/
 
   /**
    * Perform user-defined coarsening operations.  This member function
    * is called before the other coarsening operators.  The preprocess
    * function should copy data from the source components of the fine
    * patch into the source components of the destination patch on the
    * specified coarse box region.
    */
   virtual void preprocessCoarsen(
      Patch<NDIM>& coarse,
      const Patch<NDIM>& fine,
      const Box<NDIM>& coarse_box,
      const IntVector<NDIM>& ratio)
      {
         (void) coarse;
         (void) fine;
         (void) coarse_box;
         (void) ratio;
      }

   /*!*
    * Perform user-defined coarsening operations.  This member function
    * is called after the other coarsening operators.  The postprocess
    * function should copy data from the source components of the fine
    * patch into the source components of the destination patch on the
    * specified coarse box region.
    */
   virtual void postprocessCoarsen(
      Patch<NDIM>& coarse,
      const Patch<NDIM>& fine,
      const Box<NDIM>& coarse_box,
      const IntVector<NDIM>& ratio)
      {
         (void) coarse;
         (void) fine;
         (void) coarse_box;
         (void) ratio;
      }

   /*!*
    * Return stencil width of conservative averaging operations.
    */
   IntVector<NDIM> getCoarsenOpStencilWidth() const 
   {
      return(IntVector<NDIM>(0));
   }

private:
   /*!
    * These private member functions read data from input and restart.
    * When beginning a run from a restart file, all data members are read
    * from the restart file.  If the boolean flag is true when reading
    * from input, some restart values may be overridden by those in the
    * input file.
    *
    * An assertion results if the database pointer is null.
    */
   virtual void getFromInput(tbox::Pointer<tbox::Database> db);

   /*!
    * Object name used for error/warning reporting and as a label
    * for restart database entries.
    */
   string d_object_name;

   /*!
    * Variable<NDIM> - u
    */ 
   Pointer< NodeVariable<NDIM,double> > d_u_var;
   Pointer< NodeVariable<NDIM,double> > d_exact_var;
   Pointer< NodeVariable<NDIM,double> > d_error_var;

   /*!
    * Variable<NDIM> Contexts:
    *    time n+1 - new
    *    time n   - current
    *    time n-1 - old
    */ 
   Pointer<VariableContext>        d_new_cxt;
   Pointer<VariableContext>        d_cur_cxt;
   Pointer<VariableContext>        d_old_cxt;

   /*!
    * Patch<NDIM> Data ids - used to access data off patch
    */ 
   int d_u_new_id;
   int d_u_cur_id;
   int d_u_old_id;
   int d_soln_scr_id;
   int d_exact_cur_id;
   int d_error_cur_id;

   /*!
    * Refine schedule(s) used to communicate data on each level.
    */
   tbox::Array<Pointer<RefineSchedule<NDIM> > > d_bdry_fill_sched;
   
   /*!
    * Grid geometry 
    */
   Pointer<CartesianGridGeometry<NDIM> > d_grid_geometry;

   /*!
    * Initial values for constant a(x), b(x)
    */
   double d_initial_a;
   double d_initial_b;

   /*!
    * CFL condition used to determine timestep.
    */
   double d_cfl;

  /*!
   * CELL tagging variables 
   */	
	
Pointer<VariableContext> d_context_persistent;
double d_adaption_threshold;
int d_scalar_persistent;
double d_scalar;

/*!
    @brief Compute error estimator (for adaption or plotting).

    Computes in the box defined by @c estimate_data.
  */
  void computeAdaptionEstimate(
    pdat::CellData<NDIM,double> &estimate_data,
    const pdat::CellData<NDIM,double> &soln_cell_data) const;

   /*!
    * Boundary condition values.
    */
   double d_dirichlet_edge[2*NDIM];
   double d_dirichlet_node[2*NDIM];

};

