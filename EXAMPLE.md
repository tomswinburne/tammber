# Example usage of TAMMBER
Example submit scripts and configuration files can be found in the
`sample-input/vac-Fe` directory. As most parameters can
be left unchanged, here we walk through what should be altered to run `TAMMBER`
on systems of interest to you!

- [File locations](#1)
- [Specifying initial configuration(s)](#2)
- [Configuring simulation parameters in `LAMMPS`](#3)
- [State identification](#4)
- [Configuring the Markov model managing the sampling](#5)
  - [Crystal symmetry for diffusion tensor](#6)
  - [Parameters controlling accelerated sampling](#7)
  - [Cluster Monitoring](#9)
- [Analyzing output](#8)



## File locations<a name="1"></a>
When launching `TAMMBER` from some directory `execution-directory` it must
be able to find the configuration file in a directory named `input`,
i.e. the directory from where `TAMMBER` will be run **must** contain:
```bash
  execution-directory:
    input:
      ps-config.xml
```


We also need to provide an initial input configuration and interatomic potential,
which for this example are also in the `input` folder, giving for the present
example:
```bash
  execution-directory:
    input:
      ps-config.xml
      vac-Fe.dat # initial configuration
      Fe.eam.fs # interatomic potential (Marinica '07)
    sub.slurm # submit script
```

After a simulation, `TAMMBER` will write restart files (see [below](#8)) and so can be restarted from the same directory to continue sampling


## Specifying initial configuration(s)<a name="2"></a>
It is possible to seed `TAMMBER` with multiple initial configurations, but for
the moment we use only one, as specified above:
```xml
<!-- list of initial configurations; here we have one -->
<InitialConfigurations>
input/initial.dat
</InitialConfigurations>
```
## Configuring simulation parameters in `LAMMPS`<a name="3"></a>
Atomistic simulations are controlled with snippets of `LAMMPS` scripts,
so if you have used `LAMMPS` you should be able to use `TAMMBER` !

All terms in percent signs such as `%Temperature%` are fields which are filled
at runtime. Many are specified in the `TASK_XXX` tags at the bottom
of `input/ps-config.xml`

```xml
<LAMMPSEngine>
  <!-- Print logs for each LAMMPS worker. For debugging / test purposes only -->
  <LogLammps>0</LogLammps>
  <!-- Generic input script to load a configuration file %Filename%-->
  <InitScript>
    clear
    units metal
    atom_style atomic
    atom_modify map array sort 0 0.0
    dimension 3
    boundary p p p
    read_data %Filename%
    pair_style eam/fs
    pair_coeff * * input/Fe.eam.fs Fe
    neighbor 2.0 bin
    neigh_modify every 1 delay 0 check yes
    run 0
  </InitScript>

  <!-- Apply any external fixes after loading in config-->
  <PostInitScript>
  </PostInitScript>

  <!-- Run MD -->
  <MDScript>
    timestep %Timestep%

    fix NVE all nve
    fix T all langevin %Temperature% %Temperature% %Friction% %RANDU% gjf no

    thermo 10
    fix fixcom all recenter INIT INIT INIT
    run %Nsteps%
    unfix fixcom

    unfix NVE
    unfix T

  </MDScript>

  <!-- Run Minimize -->
  <MinScript>
    min_style fire
    minimize 0.0  %Tolerance% 400 400
    min_style cg
    minimize 0.0 %Tolerance% 1000 1000
  </MinScript>

  <!-- Write File -->
  <WriteScript>
    write_data %Filename%
  </WriteScript>

  <!-- Initialize Velocities -->
  <VelocityInitScript>
    velocity    all create %InitTemperature%  %RANDU% dist gaussian
  </VelocityInitScript>

</LAMMPSEngine>
```

## State identification<a name="4"></a>
`TAMMBER` identifies states by first minimizing then creating a graph of bonds
between atoms that are closer than some cutoff. For each species pair, we
thus specify a cutoff. For pure Fe, we use the `1/2<111>` bond length:
```xml
<!-- Criteria for connectivity graph -->
<Type> ConnectivityGraphStateLabeler </Type>
<Bonds>
    <Bond>
        <Between> 1 1 </Between>
        <Cutoff> 2.664 </Cutoff>
    </Bond>
</Bonds>
```


## Configuring the Markov model managing the sampling<a name="5"></a>

### Crystal symmetry for diffusion tensor<a name="6"></a>
Whilst not used during sampling, we must specify some structural properties to
allow the construction of the diffusion tensor in post processing. In particular,
we will need the lattice constant, primitive unit cell, the conventional unit cell
and the super cell.

For the primitive unit cell we can just specify the lattice type for bcc or fcc:
```xml
<!-- The following are equivalent : -->
<Lattice> 2.8552 bcc </Lattice> <!-- currently bcc or fcc only -->
```
Alternatively we can give the full cell and lattice constant as:
```xml
<Lattice>2.8552</Lattice>
<PrimitiveUnitCell>
-0.5 +0.5 +0.5
+0.5 -0.5 +0.5
+0.5 +0.5 -0.5
</PrimitiveUnitCell>
```

We then give the conventional unit cell in the same Cartesian basis but in lattice units:
```xml
<!-- Cartesian UnitCell vectors in lattice units-->
<UnitCell>
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
</UnitCell>
```


Finally, we specify the super cell vectors `[A,B,C]` from UnitCell or PrimitiveUnitCell
With UnitCell = `[a,b,c]`, give vector `v` such that `[A,B,C]=[v[0]xa,v[1]xb,v[2]xc]`
```xml
<SuperCell>7 7 7</SuperCell> <!-- i.e. A=a[700], B=a[070], C=a[007] -->
```
Alternatively, with PrimitiveUnitCell matrix `P`, give matrix `M` such that `[A,B,C]=MxP`
```xml
<SuperCell>
+0.0 +7.0 +7.0
+7.0 +0.0 +7.0
+7.0 +7.0 +0.0
</SuperCell>
```

### Parameters controlling accelerated sampling<a name="7"></a>
Perhaps the most important parameter to specify is the max temperature
for sampling. Using the available simulation data `TAMMBER` will estimate
the most suitable temperature to execute sampling within the range specified by
`<MinTemperature>` and `<MaxTemperature>`, but this currently relies on harmonic
transition state theory being valid, which in turn means no phase transitions etc!

As a result, it is best to limit `<MaxTemperature>` to approx. 80% of the melting temperature
or whatever experience with a particular system dictates, meaning for Fe we have
```xml
<!-- Lowest target temperature, in K -->
<TargetTemperature>300.0</TargetTemperature>
<!-- Temperature range for TAD sampling, steps inclusive -->
<MinTemperature>300.0</MinTemperature>
<MaxTemperature>600.0</MaxTemperature>
<TemperatureSteps>7</TemperatureSteps>
```

We also have to control how jobs are allocated.

Choice of initial distribution:
0: Delta function on first state in InitialConfigurations
1: Boltzmann on all discovered states
2: Uniform on all states (good to drive discovery for single funnel landscapes)
```xml
<RhoInitFlavor> 1 </RhoInitFlavor> <!-- i.e. Boltzmann -->
```

Scheme used to allocate workers:
0: Such that validity time increases as fast as possible (within available info)
1: Proportional to probability of finding new event from state
```xml  
<AllocScheme> 0 </AllocScheme>
```

Finally, two extra parameters which can help finish a sampling task.

As we specify a fixed simulation time, there may be incomplete `NEB` calculations
at the end of a run, meaning the constructed model is incomplete. There is thus
an option to exeute only  `TAMMBER` j
```xml
  <!--
  0: Normal operation, MD and NEB
  1: Only complete remaining NEBs
  -->
  <OnlyNEBS> 0 </OnlyNEBS>
```

### Cluster Monitoring<a name="9"></a>
For diffusion problems we typically want only one migrating object, but at high temperatures
clusters can break apart. We want to capture this breakup, but only starting sampling runs 
with one cluster. 
We use a cluster identification method to count the number of clusters, "carving" out
clusters by selecting atoms with a centrosymmetry higher than some threshold.
In the `<MarkovModel>` tags we can restrict sampling to a single cluster with
```xml
  <!-- Only sample in states with ClusterThresh clusters or less. Disabled if ClusterThresh<=0 -->
  <ClusterThresh>1</ClusterThresh>
```
We also need to specify the parameters of `TASK_CARVE` below:

```xml
<!-- CentroSymmetry based carving of clusters -->
<TaskParameter>
  <Task> TASK_CARVE </Task>
  <Flavor> 0 </Flavor> <!-- Default -->
  <!-- bcc -->
  <CentroNeighbors> 8 </CentroNeighbors>
  <!-- NB: Much higher for surfaces! -->
  <Threshold>5.0</Threshold>
  <!-- Typically ratio of 2nd neighbor length to 1st -->
  <RelativeCutoff>1.5</RelativeCutoff>
</TaskParameter>
```
The `<RelativeCutoff>` tag is a little hack to ensure the connectivity of the
"carved" out configuration can be identified using the same `<Bonds>` as
above. For our example, we know we will leave a "cage" of atoms around
the vacancy in bcc are separated by <100> (2nd nn bond length)
which is longer than the target 1/2<111> bond cutoff. We therefore scale by
|100|/|1/2(111)|~1.5 to ensure this is counted as one cluster.
Generally, RelativeCutoff ~ (2nd nn bond length) / (1st nn bond length).

## Analyzing output<a name="8"></a>
After a simulation run, you will see the following files:
```bash
  execution-directory:
    input:
      ps-config.xml
      vac-Fe.dat # initial configuration
      Fe.eam.fs # interatomic potential (Marinica '07)
    db0:
      min.db # database of configurations
      min.offsets # their location
    TammberModel.chk # restart file
    output.slurm # slurm output, if using slurm...
    err.slurm # slurm output, if using slurm...
    sub.slurm
```

To construct the input file for our python scripts, run `tammber-analyze`:
```bash
cd execution-directory
/path/to/tammber/build/tammber-analyze
```
This will produce a file `MarkovModel.xml` which is read by the `tammberAnalysis`
python package provided `process/tammberAnalysis.py`, with example usage in
`process/Diffusion_Model_example.ipynb`

To extract all configurations for your own analysis, run `tammber-dbextract`
on two cores:
```bash
cd execution-directory
mpirun -np 2 /path/to/tammber/build/tammber-dbextract
```
