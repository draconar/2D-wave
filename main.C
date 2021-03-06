//
// File:        main.C
// Package:     SAMRAI application
// Copyright:   (c) 1997-2003 The Regents of the University of California
// Release:     $Name:  $
// Revision:    $Revision: 1.7 $
// Modified:    $Date: 2003/04/01 22:06:00 $
// Description: Main program for SAMRAI wave eqn ex. problem.
//

#include "SAMRAI_config.h"

// Headers for basic SAMRAI objects
#include "tbox/SAMRAIManager.h"
#include "tbox/Database.h"
#include "tbox/InputDatabase.h"
#include "tbox/InputManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Pointer.h"
#include "tbox/PIO.h"
#include "tbox/Utilities.h"
#include "tbox/IOStream.h"
#include "VariableDatabase.h"
#include "VisItDataWriter.h"

// Headers for major algorithm/data structure objects
#include "BergerRigoutsos.h"
#include "CartesianGridGeometry.h"
#include "GriddingAlgorithm.h"
#include "StandardTagAndInitialize.h"
#include "PatchHierarchy.h"
#include "LoadBalancer.h"

// Header for application-specific algorithm/data structure object
#include "LinWaveEqn.h"
#include <string.h>
#include <stdio.h>
using namespace std;
using namespace SAMRAI;
using namespace tbox;
using namespace mesh;
using namespace geom;
using namespace appu;
using namespace hier;

/**
 * This is the main program for an example case that uses SAMRAI
 * classes to solve the 2nd order linear wave eqn (shown below in 1d):
 *
 *     u^(n+1)_i = dt^2/h^2*(u^n_i+1 - 2u^n_i + u^n_i-1) 
 *                 + 2*dt^2*u^n_i - dt^2*u^(n-1)_i
 * 
 * The main program constructs the various SAMRAI gridding objects 
 * and performs the time-stepping loop.  The program should be 
 * executed as:
 *
 *    executable <input file name>
 *
 * The input parameters for this class include:
 *
 *    max_steps       - maximum number of timesteps
 *    dt              - size of timestep
 *    end_time        - time to stop (if max_steps hasn't been reached)
 *    regrid_interval - timestep interval between regrid steps.
 *
 *  IO options:
 *
 *    log_file_name     - where logged information is written
 *    visit_dump_interval - timestep interval for writing viz dumps
 *    visit_dump_dirname  - directory name for viz dumps
 *    visit_number_procs_per_file - number of procs writing to each file.
 */

int main( int argc, char *argv[] )
{
   /*
    * Initialize MPI, SAMRAI, and enable logging.
    */

   tbox::SAMRAI_MPI::init(&argc, &argv);
   SAMRAIManager::startup();

   string input_filename;

   if ( argc != 2 ) {
        tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
           << "[options]\n" << endl;
        tbox::SAMRAI_MPI::abort();
        return (-1);
   } else {
      input_filename = argv[1];
   }

   tbox::plog << "input_filename = " << input_filename << endl;

   /****************************************************************
    *                                                              *
    *  PROBLEM SETUP                                               *
    *                                                              *
    ****************************************************************
    *                                                              *
    *  Read data from input file and initialize SAMRAI classes     *
    *                                                              *
    ****************************************************************/

   /*
    * Create input database and parse all data in input file. 
    */

   Pointer<Database> input_db = new tbox::InputDatabase("input_db");
   InputManager::getManager()->parseInputFile(input_filename,input_db);

   /*
    * Retrieve "Main" section of the input database. 
    */

   Pointer<Database> main_db = input_db->getDatabase("Main");

   int max_steps = main_db->getIntegerWithDefault("max_steps",10);
   double end_time = main_db->getDoubleWithDefault("end_time",1.);
   int regrid_interval = main_db->getIntegerWithDefault("regrid_interval",1);
   
   string log_file_name = "linwaveqn.log";
   if (main_db->keyExists("log_file_name")) {
      log_file_name = main_db->getString("log_file_name");
   }
   bool log_all_nodes = false;
   if (main_db->keyExists("log_all_nodes")) {
      log_all_nodes = main_db->getBool("log_all_nodes");
   }
   if (log_all_nodes) {
      PIO::logAllNodes(log_file_name);
   } else {
      PIO::logOnlyNodeZero(log_file_name);
   }

   int visit_dump_interval = 0;
   if (main_db->keyExists("visit_dump_interval")) {
      visit_dump_interval = main_db->getInteger("visit_dump_interval");
   } 
   string visit_dump_dirname;
   int visit_number_procs_per_file = 1;
   if ( visit_dump_interval > 0 ) {
      if (main_db->keyExists("visit_dump_dirname")) {
         visit_dump_dirname = main_db->getString("visit_dump_dirname");
      }
      if (main_db->keyExists("visit_number_procs_per_file")) {
         visit_number_procs_per_file = 
            main_db->getInteger("visit_number_procs_per_file");
      }
   }

   /*
    * The grid geometry defines the grid type (e.g. cartesian, spherical,
    * etc.).  Because SAMRAI operates on block structured indices, it can
    * support any grid geometry that may be represented as an orthogonal
    * grid.
    */
   Pointer<CartesianGridGeometry<NDIM> > grid_geometry =
      new CartesianGridGeometry<NDIM>("CartesianGeometry",
                                input_db->getDatabase("CartesianGeometry"), 1);

   /*
    * The patch hierarchy defines the adaptive grid system.
    */
   Pointer<PatchHierarchy<NDIM> > patch_hierarchy = 
      new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);

   /*
    * Set up Visualization writer.  
    */
   Pointer<VisItDataWriter<NDIM> > visit_data_writer =
      new VisItDataWriter<NDIM>("LinWaveEqn VisIt Writer",
                          visit_dump_dirname,
                          visit_number_procs_per_file);

   /*
    * This is our problem class.  See the class header for comments on it.
    */
   LinWaveEqn* wave_eqn_model = new LinWaveEqn(
      "LinWaveEqn",
      input_db->getDatabase("LinWaveEqn"),
      grid_geometry);

   /*
    * The StandardTagAndInitialize<NDIM> class performs a variety of operations
    * with user-specified parameters related to adptive gridding.  For example,
    * it manages initialization of a level, cell tagging using a gradient
    * detector, and methods to reset data after the hierarchy has been 
    * regridded.  
    *
    * It may seem overkill to use this class for this simple problem, but 
    * it has been added to facilitate adaptive gridding later.
    */
   Pointer<StandardTagAndInitialize<NDIM> > tag_and_init_ops =
      new StandardTagAndInitialize<NDIM>(
             "StandardTagAndInitialize",
             wave_eqn_model,
             input_db->getDatabase("StandardTagAndInitialize"));

   /*
    * The gridding algorithm manages adaptive gridding.  It expects a 
    * clustering scheme (i.e. how to cluster tagged-cells into patches),
    * and a load balance scheme to distribute work to processors.  In general
    * the baseline classes provided in SAMRAI should suffice for most
    * problems. It also requires a class that defines the particular tag 
    * and initialization ops that correlate with the users problem.  For
    * this, we use the "tag_and_init_ops" above, which references our 
    * "wave_eqn_model" problem class to define the user-specific operations.
    */
   Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();

   Pointer<LoadBalancer<NDIM> > load_balancer = 
      new LoadBalancer<NDIM>("LoadBalancer", input_db->getDatabase("LoadBalancer"));

   Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm = 
      new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                            input_db->getDatabase("GriddingAlgorithm"),
                            tag_and_init_ops,
                            box_generator,
                            load_balancer);

   /*
    * After creating all objects and initializing their state, we
    * print the input database and variable database contents to
    * the log file.
    */

   wave_eqn_model->setupVariables(visit_data_writer);
   tbox::plog << "\nCheck input data and variables before simulation:" << endl;
   tbox::plog << "Input database..." << endl;
   input_db->printClassData(plog);
   tbox::plog << "\nVariable database..." << endl;
   VariableDatabase<NDIM>::getDatabase()->printClassData(plog);

   

   /****************************************************************
    *                                                              *
    *  INITIALIZE DATA ON PATCHES                                  *
    *                                                              *
    ****************************************************************
    *                                                              *
    *  Build patch hierarchy and initialize the data on the patches*
    *  in the hierarchy. Note: this step is performed by the       *
    *  TimeRefinementIntegrator<NDIM> in Euler/LinAdv example cases.     *
    *  1) Create a "tag_buffer" for each level in the Hierarchy.   *
    *  2) Create the coarse (i.e. level 0) grid.                   *
    *  3) Cycle through levels 1-max_levels, initializing data     *
    *     on each.  The makeFinerLevel method calls the error      *
    *     estimator (remember, it was registered with the          *
    *     gridding algorithm object) and tags cells for refinement *
    *     as it generates patches on the finer levels.             *
    *  4) Setup communication.  Note: the communication operations *
    *     currently are very simple and not appropriately set up   *
    *     for re-gridding.  They will have to be modified for more *
    *     complex cases.                                           *
    *                                                              *
    ****************************************************************/


   tbox::Array<int> tag_buffer_array(gridding_algorithm->getMaxLevels());
   for (int il = 0; il < gridding_algorithm->getMaxLevels(); il++) {
      tag_buffer_array[il] = 1;
      tbox::pout << "il = " << il << " tag_buffer = " << tag_buffer_array[il] 
           << endl;
   }

   double loop_time = 0.;

   gridding_algorithm->makeCoarsestLevel(patch_hierarchy,loop_time);
   
   bool done = false;
   bool initial_time = true;
   for (int ln = 0; 
        gridding_algorithm->levelCanBeRefined(ln) && !done; 
        ln++) {
      gridding_algorithm->makeFinerLevel(patch_hierarchy,
                                         loop_time,
                                         initial_time,
                                         tag_buffer_array[ln]);
      done = !(patch_hierarchy->finerLevelExists(ln));
   }

   /*
    * Now that the grid hierarchy is built, setup communication for the 
    * main time advance loop.
    */
   wave_eqn_model->setupCommunication(patch_hierarchy);
   

   /*******************************************************************
    *                                                                 *
    *  MAIN TIME ADVANCE LOOP                                         *
    *                                                                 *
    *******************************************************************
    *                                                                 *
    *  1) Set start and end time.                                     *
    *  2) Start integration timesteps.                                *
    *     While (loop_time < end_time) {                              *
    *        2a) Write restart and visit data.                        *
    *        2b) Advance all levels in the hierarchy by time          *
    *            dt by calling MethodOfLinesIntegrator<NDIM>'s              *
    *            advanceHierarchy method.                             *
    *        2c) Check if it is time to do a regrid step.  If so,     *
    *            have the GriddingAlgorithm<NDIM> call its regridAllFiner   *
    *            Levels method, passing in the coarsest (0) level.    *
    *            This method will invoke the Gradient detector,       *
    *            Berger Rigoutsos algorithm, and load balancer        *
    *            while generating finer grid levels.                  *
    *     }                                                           *
    *                                                                 *
    ******************************************************************/

   wave_eqn_model->compareToExact(patch_hierarchy, loop_time);

   int iteration_num = 0;
   double dt = wave_eqn_model->getDt(patch_hierarchy);
   
   visit_data_writer->writePlotData(patch_hierarchy,
                                    iteration_num,
                                    loop_time);

   while ( (loop_time < end_time) && 
           (iteration_num < max_steps) ) {

      iteration_num++;

      tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
      tbox::pout << "At begining of timestep # " << iteration_num - 1 << endl;
      tbox::pout << "Simulation time is " << loop_time << endl;
      tbox::pout << "\n++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;

      wave_eqn_model->advance(patch_hierarchy, dt, loop_time);
 
      loop_time += dt;

      wave_eqn_model->compareToExact(patch_hierarchy, loop_time);


      tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
      tbox::pout << "At end of timestep # " << iteration_num - 1 << endl;
      tbox::pout << "Simulation time is " << loop_time << endl;
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;

      /*
       * At specified intervals, write out data files for plotting.
       * The visit_data_writer dumps data in a format that can
       * be processed by the VisIt tool. 
       */
      if ( (iteration_num % visit_dump_interval) == 0 ) {
         visit_data_writer->writePlotData(patch_hierarchy,
                                          iteration_num,
                                          loop_time);
      }

      /* 
       *  At specified intervals, regrid.
       */
      if ((iteration_num % regrid_interval) == 0 && 
           gridding_algorithm->getMaxLevels() > 1 ) {
         tbox::pout << "\n\n############################################" << endl;
         tbox::pout << "                 REGRIDDING" << endl;
         tbox::pout << "Finest level before regrid: " 
              << patch_hierarchy->getFinestLevelNumber() << endl;

         gridding_algorithm->regridAllFinerLevels(patch_hierarchy,
                                                  0,
                                                  loop_time,
                                                  tag_buffer_array);

         tbox::pout << "Finest level after regrid: " 
              << patch_hierarchy->getFinestLevelNumber() << endl;
         tbox::pout << "############################################\n\n" << endl;
      }

   }

   /*
    * At conclusion of simulation, deallocate objects.
    */

   patch_hierarchy.setNull();
   grid_geometry.setNull();

   box_generator.setNull(); 
   load_balancer.setNull();
   tag_and_init_ops.setNull();
   gridding_algorithm.setNull();
   visit_data_writer.setNull();

   //if (wave_eqn_model) delete wave_eqn_model;

   SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();
 
   return(0); 
}



