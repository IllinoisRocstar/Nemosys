# this class keeps track of definition of the mfem problem
# Al 6061t6 props taken from
# http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=ma6061t6
# Consistent units: gr mm ms MPa N N.mm
import sys
from mfem.common.arg_parser import ArgParser
from os.path import expanduser, join
import numpy as np
from numpy import sqrt, pi, cos, sin, hypot, arctan2
from scipy.special import erfc
from mfem import path
from mfem.par import intArray, add_vector
import mfem.par as mfem
from mpi4py import MPI
from vtk import vtkPoints
from vtk import vtkIdList
from vtk import vtkIntArray
from vtk import vtkDoubleArray
from vtk import vtkUnstructuredGrid
from vtk import VTK_TRIANGLE
from vtk import VTK_TETRA
from vtk import VTK_QUAD
from vtk import VTK_HEXAHEDRON
from vtk import vtkXMLUnstructuredGridWriter

# Linear elastic dynamic problem
class StiffnessOperator(mfem.PyTimeDependentOperator):
    def __init__(self, fespace,
            lmbda = 1., mu = 1., rho = 1., visc = 0.0,
            vess_tdof_list = None, vess_bdr = None,
            xess_tdof_list = None, xess_bdr = None,
            v_gfBdr = None, x_gfBdr = None,
            deform = None, velo = None,
            vx = None):
        mfem.PyTimeDependentOperator.__init__(self,
                2*fespace.GetTrueVSize(), 0.0)
        self.lmbda = lmbda
        self.mu = mu
        self.viscosity = visc
        self.deform = deform
        self.velo = velo
        self.x_gfBdr = x_gfBdr
        self.v_gfBdr = v_gfBdr
        self.vx = vx
        self.z = mfem.Vector(self.Height()/2)
        self.z.Assign(0.0)
        self.w = mfem.Vector(self.Height()/2)
        self.w.Assign(0.0)
        self.tmpVec = mfem.Vector(self.Height()/2)
        self.tmpVec.Assign(0.0)
        self.fespace = fespace
        self.xess_bdr = xess_bdr
        self.vess_bdr = vess_bdr
        self.xess_tdof_list = xess_tdof_list
        self.vess_tdof_list = vess_tdof_list

        # setting up linear form
        cv = mfem.Vector(3)
        cv.Assign(0.0)
        #self.zero_coef = mfem.ConstantCoefficient(0.0)
        self.zero_coef = mfem.VectorConstantCoefficient(cv)
        self.bx = mfem.LinearForm(self.fespace)
        self.bx.AddDomainIntegrator(
                mfem.VectorBoundaryLFIntegrator(self.zero_coef) )
        self.bx.Assemble()
        self.bv = mfem.LinearForm(self.fespace)
        self.bv.AddDomainIntegrator(
                mfem.VectorBoundaryLFIntegrator(self.zero_coef) )
        self.bv.Assemble()
        self.Bx = mfem.Vector()
        self.Bv = mfem.Vector()

        # setting up bilinear forms
        self.M = mfem.ParBilinearForm(self.fespace)
        self.K = mfem.ParBilinearForm(self.fespace)
        self.S = mfem.ParBilinearForm(self.fespace)

        self.ro = mfem.ConstantCoefficient(rho)
        self.M.AddDomainIntegrator(mfem.VectorMassIntegrator(self.ro))
        self.M.Assemble(0)
        self.M.EliminateEssentialBC(self.vess_bdr)
        self.M.Finalize(0)
        self.Mmat = self.M.ParallelAssemble()

        self.M_solver = mfem.CGSolver(
                self.fespace.GetComm())
        self.M_solver.iterative_mode = False
        self.M_solver.SetRelTol(1e-8)
        self.M_solver.SetAbsTol(0.0);
        self.M_solver.SetMaxIter(30);
        self.M_solver.SetPrintLevel(0);
        self.M_prec = mfem.HypreSmoother()
        self.M_prec.SetType(mfem.HypreSmoother.Jacobi)
        self.M_solver.SetPreconditioner(self.M_prec)
        self.M_solver.SetOperator(self.Mmat)

        lambVec = mfem.Vector(self.fespace.GetMesh().attributes.Max())
        print('Number of volume attributes : ' +
                str(self.fespace.GetMesh().attributes.Max()))
        lambVec.Assign(lmbda)
        lambVec[0] = lambVec[1]*1.0
        lambda_func = mfem.PWConstCoefficient(lambVec)
        muVec = mfem.Vector(self.fespace.GetMesh().attributes.Max())
        muVec.Assign(mu);
        muVec[0] = muVec[1]*1.0
        mu_func = mfem.PWConstCoefficient(muVec)
        self.K.AddDomainIntegrator(
            mfem.ElasticityIntegrator(lambda_func, mu_func))
        self.K.Assemble(0)
        # to set essential BC to zero value uncomment
        #self.K.EliminateEssentialBC(self.xess_bdr)
        #self.K.Finalize(0)
        #self.Kmat = self.K.ParallelAssemble()
        # to set essential BC to non-zero uncomment
        self.Kmat = mfem.HypreParMatrix()

        visc_coeff = mfem.ConstantCoefficient(visc)
        self.S.AddDomainIntegrator(
                mfem.VectorDiffusionIntegrator(visc_coeff))
        self.S.Assemble(0)
        #self.S.EliminateEssentialBC(self.vess_bdr)
        #self.S.Finalize(0)
        #self.Smat = self.S.ParallelAssemble()
        self.Smat = mfem.HypreParMatrix()

        # VX solver for implicit time-stepping
        self.VX_solver = mfem.CGSolver(
                self.fespace.GetComm())
        self.VX_solver.iterative_mode = False
        self.VX_solver.SetRelTol(1e-8)
        self.VX_solver.SetAbsTol(0.0);
        self.VX_solver.SetMaxIter(30);
        self.VX_solver.SetPrintLevel(0);
        self.VX_prec = mfem.HypreSmoother()
        self.VX_prec.SetType(mfem.HypreSmoother.Jacobi)
        self.VX_solver.SetPreconditioner(self.VX_prec)
        # setting up operators
        empty_tdof_list = intArray()
        self.S.FormLinearSystem(empty_tdof_list, self.v_gfBdr,
              self.bv, self.Smat, self.vx.GetBlock(0), self.Bv, 1)
        self.K.FormLinearSystem(empty_tdof_list, self.x_gfBdr,
              self.bx, self.Kmat, self.vx.GetBlock(1), self.Bx, 1)

    # update mass matrix
    def updMassMatMS(self, elmLst, elmCompDiagMassMat,
            rho = 1.0):
        self.M = mfem.ParBilinearForm(self.fespace)
        self.ro = mfem.ConstantCoefficient(rho)
        self.M.AddDomainIntegrator(
                mfem.VectorMassIntegrator(self.ro))
        print("in updateMass")
        self.M.Assemble(0, elmLst, elmCompDiagMassMat, True, True)
        self.M.EliminateEssentialBC(self.vess_bdr)
        self.M.Finalize(0)
        self.Mmat = self.M.ParallelAssemble()
        self.M_solver.SetOperator(self.Mmat)

    # update stiffness matrix
    def updStiffMatMS(self, elmLst, elmCompStiffMat):
        self.K = mfem.ParBilinearForm(self.fespace)
        lambVec = mfem.Vector(
                self.fespace.GetMesh().attributes.Max())
        lambVec.Assign(self.lmbda)
        lambVec[0] = lambVec[1]
        lambda_func = mfem.PWConstCoefficient(lambVec)
        muVec = mfem.Vector(
                self.fespace.GetMesh().attributes.Max())
        muVec.Assign(self.mu);
        muVec[0] = muVec[1]
        mu_func = mfem.PWConstCoefficient(muVec)
        self.K.AddDomainIntegrator(
            mfem.ElasticityIntegrator(lambda_func, mu_func))
        self.K.Assemble(0, elmLst, elmCompStiffMat, False, True)
        self.Kmat = mfem.HypreParMatrix()
        empty_tdof_list = intArray()
        self.K.FormLinearSystem(empty_tdof_list, self.x_gfBdr,
            self.bx, self.Kmat, self.vx.GetBlock(1), self.Bx, 1)

    def RecoverSln(self, vgf, xgf):
        for idx in range(self.vx.GetBlock(0).Size()):
          self.vx.GetBlock(0)[idx] = vgf[idx]
        for idx in range(self.vx.GetBlock(1).Size()):
          self.vx.GetBlock(1)[idx] = xgf[idx]

    def SetBCs(self):
        #print('setting BC')
        self.v_gfBdr.ProjectCoefficient(self.velo)
        self.x_gfBdr.ProjectCoefficient(self.deform)
        # setting up velocity BCs
        dofList = self.vess_tdof_list
        for idx in range(dofList.Size()):
          dofIdx = dofList[idx]
          self.vx.GetBlock(0)[dofIdx] = self.v_gfBdr[dofIdx]
        # setting up displacement BCs
        dofList = self.xess_tdof_list
        for idx in range(dofList.Size()):
          dofIdx = dofList[idx]
          self.vx.GetBlock(1)[dofIdx] = self.x_gfBdr[dofIdx]
          #print("Dof #"str(dofIdx) +"-> displacement =" \
          #        + str(self.vx.GetBlock(1)[dofIdx]))

    def rcvrSln(self):
        self.S.RecoverFEMSolution(self.vx.GetBlock(0), self.bv, self.v_gf)
        self.K.RecoverFEMSolution(self.vx.GetBlock(1), self.bx, self.x_gf)

    def Mult(self, vx, dvx_dt):
        sc = self.Height()/2
        v = mfem.Vector(vx, 0,  sc)
        x = mfem.Vector(vx, sc,  sc)
        dv_dt = mfem.Vector(dvx_dt, 0, sc)
        dx_dt = mfem.Vector(dvx_dt, sc,  sc)
        self.K.Mult(x, self.z)
        if (self.viscosity != 0.0):
          self.S.TrueAddMult(v,self.z)
        self.z.Neg()
        self.M_solver.Mult(self.z, dv_dt)
        dx_dt = v

    def ImplicitSolve(self, dt, vx, dvx_dt):
        sc = self.Height()/2
        v = mfem.Vector(vx, 0,  sc)
        x = mfem.Vector(vx, sc,  sc)
        dv_dt = mfem.Vector(dvx_dt, 0, sc)
        dx_dt = mfem.Vector(dvx_dt, sc,  sc)
        #print("v = ")
        #v.Print()
        #print("x = ")
        #x.Print()
        # By eliminating kx from the coupled system:
        # kv = -M^{-1}*[H(x + dt*kx) + S*(v + dt*kv)]
        # kx = v + dt*kv
        # we reduce it to a linear equation for kv,
        # represented by the backward_euler_oper
        ndt = -1*dt
        self.VX = mfem.Add(1.0, self.Mmat,
                ndt, self.Smat)
        self.VX_solver.SetOperator(self.VX)

        #add_vector(x, dt, dv_dt, self.w)
        self.Kmat.Mult(x, self.z)
        #self.S.TrueAddMult(v, self.z)
        self.Smat.Mult(v, self.tmpVec)
        add_vector(self.z, 1.0, self.tmpVec, self.z)
        self.z.Neg()
        self.z += self.Bx
        self.z += self.Bv
        self.VX_solver.Mult(self.z, dv_dt)
        add_vector(v, dt, dv_dt, dx_dt)
        #print("dv_dt =")
        #dv_dt.Print()
        #print("dx_dt =")
        #dx_dt.Print()

class defBCs(mfem.VectorPyCoefficientT):
   def EvalValue(self, x, t):
       dim = len(x)
       disp = np.zeros(len(x))
       if (x[0] >= 8.0):
         if (t <= 5e12):
           disp[1] = 0.0001*t
         else:
           disp[1] = 0.0
       return disp

class velBCs(mfem.VectorPyCoefficientT):
   def EvalValue(self, x, t):
       dim = len(x)
       vel = np.zeros(len(x))
       if (x[0]>= 8.0 and t<=5e12):
         vel[1]  = 0.0001
       return vel

##################################################
# orchestrator class : performs all sim managements
##################################################

class mfemSim(object):
    def __init__(self, desc = ''):
      # 1. Object initialization
      self.desc =  desc
      self.myId = MPI.COMM_WORLD.rank
      self.nProc = MPI.COMM_WORLD.size

    def __del__(self):
      # do nothing for now
      self.finalize()
      print('MFEM problem is deleted')

    # problem initialization
    def initialize(self, inMeshObj=None, inMeshFile=None):
      # 2. Problem initialization
      self.parser = ArgParser(description='Based on MFEM Ex16p')
      self.parser.add_argument('-m', '--mesh',
            default = 'beam-tet.mesh',
            action = 'store', type = str,
            help='Mesh file to use.')
      self.parser.add_argument('-rs', '--refine-serial',
            action = 'store', default = 1, type=int,
            help = "Number of times to refine the mesh \
                    uniformly in serial")
      self.parser.add_argument('-rp', '--refine-parallel',
            action = 'store', default = 0, type=int,
            help = "Number of times to refine the mesh \
            uniformly in parallel")
      self.parser.add_argument('-o', '--order',
            action = 'store', default = 1, type=int,
            help = "Finite element order (polynomial \
                    degree)");
      self.parser.add_argument('-s', '--ode-solver',
            action = 'store', default = 3, type = int,
            help = '\n'.join(
            ["ODE solver: 1 - Backward Euler, 2 - SDIRK2, \
              3 - SDIRK3", "\t\t 11 - Forward Euler, \
              12 - RK2, 13 - RK3 SSP, 14 - RK4."]))
      self.parser.add_argument('-t', '--t-final',
            action = 'store', default = 20., type=float,
            help = "Final time; start time is 0.")
      self.parser.add_argument("-dt", "--time-step",
            action = 'store', default = 5e-3, type=float,
            help = "Time step.")
      self.parser.add_argument("-v", "--viscosity",
            action = 'store', default = 0.00, type=float,
            help = "Viscosity coefficient.")
      self.parser.add_argument('-L', '--lmbda',
            action = 'store', default = 1.e0, type=float,
            help = 'Lambda of Hooks law')
      self.parser.add_argument('-mu', '--shear-modulus',
            action = 'store', default = 1.e0, type=float,
            help = 'Shear modulus for Hooks law')
      self.parser.add_argument('-rho', '--density',
            action = 'store', default = 1.0, type=float,
            help = 'mass density')
      self.parser.add_argument('-vis', '--visualization',
            action = 'store_true',
            help='Enable GLVis visualization')
      self.parser.add_argument('-vs',
            '--visualization-steps',
            action = 'store', default = 25, type = int,
            help = "Visualize every n-th timestep.")
      args = self.parser.parse_args()
      self.ser_ref_levels = args.refine_serial
      self.par_ref_levels = args.refine_parallel
      self.order = args.order
      self.dt = args.time_step
      self.visc = args.viscosity
      self.t_final = args.t_final
      self.lmbda = args.lmbda
      self.mu = args.shear_modulus
      self.rho = args.density
      self.visualization = args.visualization
      self.ti = 1
      self.vis_steps = args.visualization_steps
      self.ode_solver_type = args.ode_solver
      self.t = 0.0;
      self.last_step = False
      if self.myId == 0: self.parser.print_options(args)

      # 3. Reading mesh
      if inMeshObj is None:
        self.meshFile = inMeshFile
        if self.meshFile is None:
          self.meshFile = args.mesh
        self.mesh = mfem.Mesh(self.meshFile, 1,1)
      else:
        self.mesh = inMeshObj
      self.dim = self.mesh.Dimension()
      print("Mesh dimension: %d" % self.dim)
      print("Number of vertices in the mesh: %d " %
              self.mesh.GetNV())
      print("Number of elements in the mesh: %d " %
              self.mesh.GetNE())

      # 4. Define the ODE solver used for time integration.
      #    Several implicit singly diagonal implicit
      #    Runge-Kutta (SDIRK) methods, as well as
      #    explicit Runge-Kutta methods are available.
      if   self.ode_solver_type == 1:
          self.ode_solver = BackwardEulerSolver()
      elif self.ode_solver_type == 2:
          self.ode_solver = mfem.SDIRK23Solver(2)
      elif self.ode_solver_type == 3:
          self.ode_solver = mfem.SDIRK33Solver()
      elif self.ode_solver_type == 11:
          self.ode_solver = ForwardEulerSolver()
      elif self.ode_solver_type == 12:
          self.ode_solver = mfem.RK2Solver(0.5);
      elif self.ode_solver_type == 13:
          self.ode_solver = mfem.RK3SSPSolver()
      elif self.ode_solver_type == 14:
          self.ode_solver = mfem.RK4Solver()
      elif self.ode_solver_type == 22:
          self.ode_solver = mfem.ImplicitMidpointSolver()
      elif self.ode_solver_type == 23:
          self.ode_solver = mfem.SDIRK23Solver()
      elif self.ode_solver_type == 24:
          self.ode_solver = mfem.SDIRK34Solver()
      else:
        print("Unknown ODE solver type: "
              + str(self.ode_solver_type))
        exit

      # 5. Refine the mesh in serial to increase the
      #    resolution. In this example we do
      #    'ser_ref_levels' of uniform refinement, where
      #    'ser_ref_levels' is a command-line parameter.
      for lev in range(self.ser_ref_levels):
        self.mesh.UniformRefinement()

      # 6. Define a parallel mesh by a partitioning of
      #    the serial mesh. Refine this mesh further
      #    in parallel to increase the resolution. Once the
      #    parallel mesh is defined, the serial mesh can
      #    be deleted.
      self.pmesh = mfem.ParMesh(MPI.COMM_WORLD, self.mesh)
      for lev in range(self.par_ref_levels):
        self.pmesh.UniformRefinement()

      # 7. Define the vector finite element space
      #    representing the current and the
      #    initial temperature, u_ref.
      self.fe_coll = mfem.H1_FECollection(
           self.order, self.dim)
      self.fespace = mfem.ParFiniteElementSpace(
           self.pmesh, self.fe_coll, self.dim)
      self.fe_size = self.fespace.GlobalTrueVSize();
      if self.myId == 0:
        print( "FE Number of unknowns: " +
                str(self.fe_size))
      true_size = self.fespace.TrueVSize()
      self.true_offset = mfem.intArray(3)
      self.true_offset[0] = 0
      self.true_offset[1] = true_size
      self.true_offset[2] = 2*true_size
      self.vx = mfem.BlockVector(self.true_offset)
      self.v_gf  = mfem.ParGridFunction(self.fespace)
      self.v_gfbnd  = mfem.ParGridFunction(self.fespace)
      self.x_gf  = mfem.ParGridFunction(self.fespace)
      self.x_gfbnd  = mfem.ParGridFunction(self.fespace)
      self.x_ref = mfem.ParGridFunction(self.fespace)
      self.pmesh.GetNodes(self.x_ref)

      # 8. Set the initial conditions for u.
      #self.velo = InitialVelocity(self.dim)
      self.velo = velBCs(self.dim)
      #self.deform =  InitialDeformation(self.dim)
      self.deform = defBCs(self.dim)
      self.v_gf.ProjectCoefficient(self.velo)
      self.v_gfbnd.ProjectCoefficient(self.velo)
      self.x_gf.ProjectCoefficient(self.deform)
      self.x_gfbnd.ProjectCoefficient(self.deform)
      #self.v_gf.GetTrueDofs(self.vx.GetBlock(0));
      #self.x_gf.GetTrueDofs(self.vx.GetBlock(1));

      # setup boundary-conditions
      self.xess_bdr = mfem.intArray(
              self.fespace.GetMesh().bdr_attributes.Max())
      self.xess_bdr.Assign(0)
      self.xess_bdr[0] = 1;
      self.xess_bdr[1] = 1;
      self.xess_tdof_list = intArray()
      self.fespace.GetEssentialTrueDofs(self.xess_bdr,
              self.xess_tdof_list)
      #print('True x essential BCs are')
      #self.xess_tdof_list.Print()

      self.vess_bdr = mfem.intArray(
              self.fespace.GetMesh().bdr_attributes.Max())
      self.vess_bdr.Assign(0)
      self.vess_bdr[0] = 1;
      self.vess_bdr[1] = 1;
      self.vess_tdof_list = intArray()
      self.fespace.GetEssentialTrueDofs(self.vess_bdr,
              self.vess_tdof_list)
      #print('True v essential BCs are')
      #self.vess_tdof_list.Print()

      # 9. Initialize the stiffness operator
      self.oper = StiffnessOperator(self.fespace,
              self.lmbda, self.mu,
              self.rho, self.visc,
              self.vess_tdof_list, self.vess_bdr,
              self.xess_tdof_list, self.xess_bdr,
              self.v_gfbnd, self.x_gfbnd,
              self.deform, self.velo,
              self.vx)

      # 10. Setting up file output
      self.smyid = '{:0>2d}'.format(self.myId)

      # initializing ode solver
      self.ode_solver.Init(self.oper)

    # stepping in time
    def step(self, nStep):
      # initialization
      self.last_step = False
      tFinal = self.t + self.dt*nStep
      if (self.t + self.dt >= tFinal - self.dt/2):
          self.last_step = True
      # loop for time stepping
      while not self.last_step:
        if (self.t + self.dt >= tFinal - self.dt/2):
            self.last_step = True
        # update BCs
        self.deform.SetTime(self.t)
        self.velo.SetTime(self.t)
        # step in time
        self.t, self.dt = self.ode_solver.Step(
                self.vx, self.t, self.dt)
        # writing to output
        if (self.last_step or (self.ti % self.vis_steps) == 0):
          if self.myId == 0:
            print("Running step " + str(self.ti)
                  + ", t = " + str(self.t))
          self.oper.SetBCs()
          self.v_gf.Distribute(self.vx.GetBlock(0))
          self.x_gf.Distribute(self.vx.GetBlock(1))
          self.solToFile()
        # preparing for the next time step
        self.ti = self.ti + 1
        self.oper.SetBCs()

    # sets up simulation agent after a major update
    # such as mass matrix or stiffness matrix
    def postUpdate(self):
      self.oper.RecoverSln(self.v_gf, self.x_gf)

    # finalize the simulation
    def finalize(self):
      print('mfem simulation is finalizing')

    # writes solution to file
    def solToFile(self):
      # print mesh to vtk format
      sol_name   =  "fe-sol-" + self.smyid \
      + '-' + str(self.ti).zfill(6) + '.vtu'
      #sol_name   =  "fe-sol-" + self.smyid \
      #+ '-' + str(self.ti).zfill(6) 
      #self.vtkStream = open(sol_name+'.vtk', 'w')
      #self.pmesh.PrintVTK(self.vtkStream, 1, 1)
      fieldName = 'V'
      #self.v_gf.SaveVTK(self.vtkStream, str(fieldName), 1)
      fieldName = 'U'
      #self.x_gf.SaveVTK(self.vtkStream, str(fieldName), 1)
      #self.vtkStream.close()
      meshToVtu(self.pmesh, sol_name, (self.v_gf,self.x_gf), ('V','U'))

    # get solution vector
    def getSolVec(self, name):
      ret = []
      name = name.lower()
      if name == 'displacement':
        sz = self.vx.GetBlock(1).Size()
        ret = np.zeros(sz, np.double)
        for i in range(sz):
          ret[i] = self.vx.GetBlock(1)[i]
      elif name == 'velocity':
        sz = self.vx.GetBlock(0).Size()
        ret = np.zeros(sz, np.double)
        for i in range(sz):
          ret[i] = self.vx.GetBlock(0)[i]
      else:
        raise ValueError('Requested solution is invalid.')
      return(ret)


def meshToVtu(mesh, out_mesh, grid_func, grid_func_names):
#  mesh.UniformRefinement()
  # get mesh dimensions (1,2,3)
  dim = mesh.Dimension()
  # read nodes from mesh and load into vtk points 
  nodes = mfem.Vector()
  mesh.GetNodes(nodes) 
  numPoints = nodes.Size()/dim
  points = vtkPoints()
  points.SetNumberOfPoints(numPoints)
  for j in range(0,numPoints):
    if (dim == 2):
      points.InsertPoint(j,nodes[j],nodes[j+numPoints],0)
    elif (dim == 3):
      points.InsertPoint(j,nodes[j],nodes[j+numPoints],nodes[j+2*numPoints])
  
  # put points in vtkunsgrid
  vtkmesh = vtkUnstructuredGrid()
  vtkmesh.SetPoints(points)
   
  # read boundary elements and attributes from mesh and load into vtk cells 
  numBndElem = mesh.GetNBE()
  numElem = mesh.GetNE()
  attributes = vtkIntArray()
  attributes.SetNumberOfComponents(1) 
  attributes.SetNumberOfTuples(numBndElem+numElem)
  attributes.SetName("material")
  for j in range(0,numBndElem):
    mfem_cell = mesh.GetBdrElementVertices(j)
    vtk_cell = vtkIdList()  
    for k in range(0,len(mfem_cell)):
      vtk_cell.InsertNextId(mfem_cell[k])
    cell_type = mesh.GetBdrElementType(j)
    if cell_type == 2:
      vtkmesh.InsertNextCell(VTK_TRIANGLE, vtk_cell) 
    elif cell_type == 3:
      vtkmesh.InsertNextCell(VTK_QUAD, vtk_cell)   
    attributes.SetTuple1(j,mesh.GetBdrAttribute(j))

  # read interior elements and attributes from mesh and load into vtk cells
  for j in range(0,numElem):
    mfem_cell = mesh.GetElementVertices(j)
    vtk_cell = vtkIdList()  
    for k in range(0,len(mfem_cell)):
      vtk_cell.InsertNextId(mfem_cell[k])
    cell_type = mesh.GetElementType(j)
    if cell_type == 4:
      vtkmesh.InsertNextCell(VTK_TETRA, vtk_cell) 
    elif cell_type == 5:
      vtkmesh.InsertNextCell(VTK_HEXAHEDRON, vtk_cell)   
    attributes.SetTuple1(j+numBndElem, mesh.GetAttribute(j))
  
  #nemMesh = meshBase.Create(out_mesh)
  #print(str(nemMesh.getNumberOfPoints()) + ' ' + str(nemMesh.getNumberOfCells()))

  if (grid_func != 0 and grid_func_names != 0):
    for func, name in zip(grid_func, grid_func_names):
      point_data = vtkDoubleArray()
      point_data.SetNumberOfComponents(func.VectorDim())
      point_data.SetName(name)
      if func.VectorDim() == 3:
        nodal_values1 = mfem.Vector()
        nodal_values2 = mfem.Vector()
        nodal_values3 = mfem.Vector()
        func.GetNodalValues(nodal_values1,1)
        func.GetNodalValues(nodal_values2,2)
        func.GetNodalValues(nodal_values3,3)
        for j in range(0,numPoints):
          point_data.InsertNextTuple3(nodal_values1[j], nodal_values2[j], nodal_values3[j])
      vtkmesh.GetPointData().AddArray(point_data)
   
  # add attributes to vtk dataSet
  vtkmesh.GetCellData().AddArray(attributes)
  writer = vtkXMLUnstructuredGridWriter()
  writer.SetFileName(out_mesh)
  writer.SetInputData(vtkmesh)
  writer.Write()
