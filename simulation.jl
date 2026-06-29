########################
#This code solves the non-linear wave equation in 2+1. 
#It implements a RK4 to evolve in time the system, 
#and fourth order finite differences for the spatial derivatives
########################
#Author: Alejandro Cardenas-Avendano

include("4dCartoonwave.jl")
include("IC.jl")
include("DOperators.jl")
include("evolution.jl")

using HDF5
using DelimitedFiles
using ArgParse
using GSL

function parse_commandline()
    s = ArgParseSettings()

    #Arguments
    @add_arg_table s begin
        "--res"
            help = "Resolution of the grid"
            arg_type = Float64
            default = 0.01
        "--ti"
            help = "Initial time"
            arg_type = Float64
            default = 0.0
        "--tf"
            help = "Final time"
            arg_type = Float64
            default = 1.0
        "--energb"
            #What percentage of the grid the energy is going to be measured
            help = "Energy bound"
            arg_type = Float64
            default = 0.95
        "--amp"
            help = "Amplitude"
            arg_type = Float64
            default = 4000.0 #Default for NLMWP, cubic nonlinearity
        "--space" #@Truong
            help = "Spacetime geometry"
            arg_type = String
            default = "NLMWP"
        "--wave" #@Truong
            help = "Wave equation"
            arg_type = String
            default = "Cubic" #cubic: phi^3
        "--quasi" #@Truong
            help = "Inclusion of quasilinear term in time derivative nonlinearity"
            arg_type = Float64
            default = 0.0 #0.0=does not include quasi term, 1.0=includes it
        "--courant" #@Truong
            help = "Courant factor"
            arg_type = Float64
            default = 0.5
        "--diss" #@Truong
            help = "Kreiss-Oliger dissipation"
            arg_type = Float64
            default = 0.3
        "--snapdt" #@Truong
            help = "Increments to store the field"
            arg_type = Float64
            default = 0.1
    end
    return parse_args(s)
end

parsed_args = parse_commandline()

#Spacetime choices: NLMWP, Minkowski, Hayward, or Bardeen @Truong
spacetime=parsed_args["space"];

#Throw an error message if an invalid spacetime is inputted @Truong
if spacetime != "NLMWP" && spacetime != "Minkowski" && spacetime != "Hayward" && spacetime != "Bardeen"
    error("Invalid spacetime geometry! Options: NLMWP, Minkowski, Hayward, Bardeen")
end

nlmwp=0.0;
hayward=0.0;
bardeen=0.0;

if spacetime == "NLMWP"
    nlmwp=1.0;
end
if spacetime == "Hayward"
    hayward=1.0;
end
if spacetime == "Bardeen"
    bardeen=1.0;
end

#Wave equation choices: Linear, Cubic (phi^3), or TimeDerivs @Truong
waveeqn=parsed_args["wave"];

#Throw an error message if an invalid wave equation is inputted @Truong
if waveeqn != "Linear" && waveeqn != "Cubic" && waveeqn != "TimeDerivs"
    error("Invalid wave equation! Options: Linear, Cubic, TimeDerivs")
end

cubic=0.0;
timederivs=0.0;

if waveeqn == "Cubic"
    cubic=1.0;
end
if waveeqn == "TimeDerivs"
    timederivs=1.0;
end

#Quasilinear term choices: 0.0 or 1.0
quasi=parsed_args["quasi"]; #@Truong

#Throw an error message if input is not 0.0 or 1.0
if quasi != 0.0 && quasi != 1.0
    error("Invalid quasilinear term! Options: 0.0 or 1.0")
end

ICtype="Analytical"

orderFD=convert(Int64,4)


#For the piecewise poly
const Am=parsed_args["amp"]; #amplitude

#Default initial data values for NLMWP
const r0=0.02;
const r1=0.6;
if spacetime=="Hayward" #@Truong
    const r1=0.66; #This is optimized for Hayward values l=0.18, m=0.20, l/m=0.9
end

const ells=[1,2] #Spherical harmonics, l modes
const polylogexp=4.0

# Resolution of the grids
const dy=parsed_args["res"];
const dz=parsed_args["res"];

# Courant factor: dt/dx
const lambda = parsed_args["courant"]; #@Truong
# Kreiss-Oliger dissipation
const epsKO = parsed_args["diss"]; #@Truong
# Initial Time
const Ti= parsed_args["ti"];
# Final Time
const Tf= parsed_args["tf"];

boundenerg=parsed_args["energb"];

#@Truong
println("\n================================")
println("Spacetime\t=\t", spacetime)
println("Wave equation\t=\t", waveeqn)
println("Amplitude\t=\t", Am)
println("Grid resolution\t=\t", dy)
println("Courant factor\t=\t", lambda)
println("Kreiss-Oliger\t=\t", epsKO)
println("Initial time\t=\t", Ti)
println("Final time\t=\t", Tf)
println("Quasilinear\t=\t", quasi)
println("Energy bound\t=\t", boundenerg)
println("Initial data r0\t=\t", r0)
println("Initial data r1\t=\t", r1)
println("================================\n")

#Spatial limits in compactified coordinates
ymin       =  0.0;
ymax       =  1.0;

zmin       =  -1.0;
zmax       =  1.0;

# Define the folder name
case_name = "$(spacetime)_$(waveeqn)_$(Am)_$(dy)_$(lambda)_$(epsKO)" #@Truong
folder_name = "Results_"*case_name #@Truong

# Check if the folder exists
if !isdir(folder_name)
    # Create the folder if it doesn't exist
    mkdir(folder_name)
    println("Folder '"*folder_name*"' created.\n") #@Truong
else
    println("Folder '"*folder_name*"' already exists.\n") #@Truong
end

location=folder_name*"/" #@Truong 

#@Truong
#For the metric potential
#Hayward
#These values allow for trapping, but the spacetime does not have a BH
# const l=0.15; #l/m=0.833
# const m=0.18;
const l=0.18; #l/m=0.9
const m=0.20;

#Bardeen
#These values allow for trapping, but the spacetime does not have a BH
const mBD=0.32;
const qBD=0.25;

#NLMWP
const a=0.026;
const b=11.20;

if spacetime=="Hayward"
    println("Using Hayward metric values: l = $(l), m = $(m)\n")
elseif spacetime=="Bardeen"
    println("Using Bardeen metric values: mBD = $(mBD), qBD = $(qBD)\n")
elseif spacetime=="Minkowski"
    const a=1e8;
    const b=1e8;
    println("Using Minkowski metric\n")
else
    println("Using NLMWP metric values: a = $(a), b = $(b)\n")
end

const dt=lambda*dy;

const Ny0= convert(Int64,(ymax-ymin)/dy+1);
const Nz0= convert(Int64,(zmax-zmin)/dz+1);


println("Allocating memory for the matrices of size $Ny0 x $Nz0");

const Nt0= convert(Int64,round((Tf-Ti)/dt));

const ts=collect(range(Ti, Tf, length=Nt0+1)); #timesteps to evolve

const roundfact=1e8

#Cadence to auxiliar quantites (energy)
const dtRaux=0.05;
#Cadence to store the field and its time derivative (to save a snapshot)
const dtRaux2=parsed_args["snapdt"]; #@Truong

const NtR=convert(Int64,round((Tf-Ti)/dtRaux));
#Cadence to save the max vals
const dtR=round(dtRaux*roundfact);
#Cadence to save images of the simulation
const dtRFS=round(dtRaux2*roundfact);

const ts_aux=round.(roundfact*ts);

#Initializations

const ys=collect(range(ymin, ymax, length=Ny0)); #0.0 to 1.0
const zs=collect(range(zmin, zmax, length=Nz0)); #-1.0 to 1.0

const pidotmaxval=collect(range(Ti, Tf, length=Nt0+1));

const phidxdxmaxval=collect(range(Ti, Tf, length=Nt0+1));
const phidydymaxval=collect(range(Ti, Tf, length=Nt0+1));
const phidzdzmaxval=collect(range(Ti, Tf, length=Nt0+1));
const phidydzmaxval=collect(range(Ti, Tf, length=Nt0+1));

#Spherical Coordinates

#Max values of second angular derivatives
const phidthdthmaxval=collect(range(Ti, Tf, length=Nt0+1)); #theta
const phidphidphimaxval=collect(range(Ti, Tf, length=Nt0+1)); #phi

yvaltest = round(Int, (boundenerg - ymin)/dy) + 1;
zvaltest1 = round(Int, (-boundenerg - zmin)/dz) + 1;
zvaltest2 = round(Int, ( boundenerg - zmin)/dz) + 1;

# yvaltest=findall(x -> x == boundenerg, ys)[1];
# zvaltest1=findall(x -> x == -boundenerg, zs)[1];
# zvaltest2=findall(x -> x == boundenerg, zs)[1];

const tsR=collect(range(Ti, Tf, length=NtR+1));

const phi_M1= zeros(Float64, (Ny0,Nz0));
const phi_M2= zeros(Float64, (Ny0,Nz0));

const phidx_M= zeros(Float64, (Ny0,Nz0));
const phidy_M= zeros(Float64, (Ny0,Nz0));
const phidz_M= zeros(Float64, (Ny0,Nz0));

const phidxdx_M= zeros(Float64, (Ny0,Nz0));
const phidxdy_M= zeros(Float64, (Ny0,Nz0));
const phidxdz_M= zeros(Float64, (Ny0,Nz0));
const phidydy_M= zeros(Float64, (Ny0,Nz0));
const phidydz_M= zeros(Float64, (Ny0,Nz0));
const phidzdz_M= zeros(Float64, (Ny0,Nz0));

const pi_M1= zeros(Float64, (Ny0,Nz0));
const pi_M2= zeros(Float64, (Ny0,Nz0));

const pidx_M= zeros(Float64, (Ny0,Nz0));
const pidy_M= zeros(Float64, (Ny0,Nz0));
const pidz_M= zeros(Float64, (Ny0,Nz0));

const phidot_M= zeros(Float64, (Ny0,Nz0));
const pidot_M= zeros(Float64, (Ny0,Nz0));

const phik1_M= zeros(Float64, (Ny0,Nz0));
const pik1_M= zeros(Float64, (Ny0,Nz0));

const phik2_M= zeros(Float64, (Ny0,Nz0));
const pik2_M= zeros(Float64, (Ny0,Nz0));

const phik3_M= zeros(Float64, (Ny0,Nz0));
const pik3_M= zeros(Float64, (Ny0,Nz0));

const phik4_M= zeros(Float64, (Ny0,Nz0));
const pik4_M= zeros(Float64, (Ny0,Nz0));

#Contravariant metric components
const gtt=zeros(Float64, (Ny0,Nz0));

const gxx=zeros(Float64, (Ny0,Nz0));
const gxy=zeros(Float64, (Ny0,Nz0));
const gxz=zeros(Float64, (Ny0,Nz0));

const gyy=zeros(Float64, (Ny0,Nz0));
const gyz=zeros(Float64, (Ny0,Nz0));

const gzz=zeros(Float64, (Ny0,Nz0));

#Derivatives of metric components

const gtt_dx=zeros(Float64, (Ny0,Nz0));
const gtt_dy=zeros(Float64, (Ny0,Nz0));
const gtt_dz=zeros(Float64, (Ny0,Nz0));

const gxx_dx=zeros(Float64, (Ny0,Nz0));
const gxx_dy=zeros(Float64, (Ny0,Nz0));
const gxx_dz=zeros(Float64, (Ny0,Nz0));

const gxy_dx=zeros(Float64, (Ny0,Nz0));
const gxy_dy=zeros(Float64, (Ny0,Nz0));
const gxy_dz=zeros(Float64, (Ny0,Nz0));

const gxz_dx=zeros(Float64, (Ny0,Nz0));
const gxz_dy=zeros(Float64, (Ny0,Nz0));
const gxz_dz=zeros(Float64, (Ny0,Nz0));

const gyy_dx=zeros(Float64, (Ny0,Nz0));
const gyy_dy=zeros(Float64, (Ny0,Nz0));
const gyy_dz=zeros(Float64, (Ny0,Nz0));

const gyz_dx=zeros(Float64, (Ny0,Nz0));
const gyz_dy=zeros(Float64, (Ny0,Nz0));
const gyz_dz=zeros(Float64, (Ny0,Nz0));

const gzz_dx=zeros(Float64, (Ny0,Nz0));
const gzz_dy=zeros(Float64, (Ny0,Nz0));
const gzz_dz=zeros(Float64, (Ny0,Nz0));

#Minus square root of the determinant of the metric

const sqrtming=zeros(Float64, (Ny0,Nz0));

const sqrtming_dx=zeros(Float64, (Ny0,Nz0));
const sqrtming_dy=zeros(Float64, (Ny0,Nz0));
const sqrtming_dz=zeros(Float64, (Ny0,Nz0));

phidotdy_M=zeros(Float64, (Ny0,Nz0))
phidotdz_M=zeros(Float64, (Ny0,Nz0))

println("Evaluating the metric and the initial conditions");
ICs!(phi_M1,pi_M1,gtt,gxx,gxy,gxz,gyy,gyz,gzz,sqrtming)

println("Computing the metric derivatives");
metricderivatives!(gtt_dx,gxx_dx,gxy_dx,gxz_dx,gyy_dx,gyz_dx,gzz_dx,sqrtming_dx,gtt_dy,gxx_dy,gxy_dy,gxz_dy,gyy_dy,gyz_dy,gzz_dy,sqrtming_dy,gtt_dz,gxx_dz,gxy_dz,gxz_dz,gyy_dz,gyz_dz,gzz_dz,sqrtming_dz)

#Stores metric components to check if needed
if isfile(location*"Metric_"*case_name*".h5")
    rm(location*"Metric_"*case_name*".h5")
end
println("Saving metric components")
h5write(location*"Metric_"*case_name*".h5", "gtt",gtt[:,:,:])
h5write(location*"Metric_"*case_name*".h5", "gxx",gxx[:,:,:])
h5write(location*"Metric_"*case_name*".h5", "gxy",gxy[:,:,:])
h5write(location*"Metric_"*case_name*".h5", "gxz",gxz[:,:,:])
h5write(location*"Metric_"*case_name*".h5", "gyy",gyy[:,:,:])
h5write(location*"Metric_"*case_name*".h5", "gyz",gyz[:,:,:])
h5write(location*"Metric_"*case_name*".h5", "gzz",gzz[:,:,:])
h5write(location*"Metric_"*case_name*".h5", "sqrtming",sqrtming[:,:,:])

#Stores first derivatives wrt x of metric components to check if needed
if isfile(location*"MetricXDerivatives_"*case_name*".h5")
    rm(location*"MetricXDerivatives_"*case_name*".h5")
end
println("Saving first derivatives wrt x of metric components")
h5write(location*"MetricXDerivatives_"*case_name*".h5", "gtt_dx",gtt_dx[:,:,:])
h5write(location*"MetricXDerivatives_"*case_name*".h5", "gxx_dx",gxx_dx[:,:,:])
h5write(location*"MetricXDerivatives_"*case_name*".h5", "gxy_dx",gxy_dx[:,:,:])
h5write(location*"MetricXDerivatives_"*case_name*".h5", "gxz_dx",gxz_dx[:,:,:])
h5write(location*"MetricXDerivatives_"*case_name*".h5", "gyy_dx",gyy_dx[:,:,:])
h5write(location*"MetricXDerivatives_"*case_name*".h5", "gyz_dx",gyz_dx[:,:,:])
h5write(location*"MetricXDerivatives_"*case_name*".h5", "gzz_dx",gzz_dx[:,:,:])
h5write(location*"MetricXDerivatives_"*case_name*".h5", "sqrtming_dx",sqrtming_dx[:,:,:])

fileInfo=location*"Info_Wave_"*case_name*".txt" #@Truong
fileMax=location*"Max_Wave_"*case_name*"_$(Ti).h5" #@Truong

if isfile(fileInfo)
    rm(fileInfo)
    println("Info File Overwritten")
end

#Save simulation parameters in Info_Wave file @Truong
if spacetime=="Hayward"
    write(fileInfo, 
"Simulation parameters
================================
Spacetime = $(spacetime)
Wave equation = $(waveeqn)
Amplitude = $(Am)
Grid resolution = $(dy)
Courant factor = $(lambda)
Kreiss-Oliger = $(epsKO)
Quasilinear = $(quasi)
Energy bound = $(boundenerg)
r0 = $(r0)
r1 = $(r1)
ells = $(ells)
polylogexp = $(polylogexp)
\n$(spacetime) metric values:
l = $(l)
m = $(m)
================================\n"
    )
elseif spacetime=="Bardeen"
    write(fileInfo, 
"Simulation parameters
================================
Spacetime = $(spacetime)
Wave equation = $(waveeqn)
Amplitude = $(Am)
Grid resolution = $(dy)
Courant factor = $(lambda)
Kreiss-Oliger = $(epsKO)
Quasilinear = $(quasi)
Energy bound = $(boundenerg)
r0 = $(r0)
r1 = $(r1)
ells = $(ells)
polylogexp = $(polylogexp)
\n$(spacetime) metric values:
mBD = $(mBD)
qBD = $(qBD)
================================\n"
    )
else
    #Minkowski, NLMWP
    write(fileInfo, 
"Simulation parameters
================================
Spacetime = $(spacetime)
Wave equation = $(waveeqn)
Amplitude = $(Am)
Grid resolution = $(dy)
Courant factor = $(lambda)
Kreiss-Oliger = $(epsKO)
Quasilinear = $(quasi)
Energy bound = $(boundenerg)
r0 = $(r0)
r1 = $(r1)
ells = $(ells)
polylogexp = $(polylogexp)
\n$(spacetime) metric values:
a = $(a)
b = $(b)
================================\n"
    )
end

println("Evolving the wave equation");
@time simulation!(size(ts)[1],dtR,phi_M1,pi_M1,phi_M2,pi_M2,phidx_M,pidx_M,phidy_M,pidy_M,phidz_M,pidz_M,phidot_M,pidot_M,phidxdx_M,phidxdy_M,phidxdz_M,phidydy_M,phidydz_M,phidzdz_M)

if isfile(fileMax)
    rm(fileMax)
    println("File Maxs Overwritten")
end
h5write(fileMax, "ts",ts[:]) # time
h5write(fileMax, "pidotmaxs",pidotmaxval[:]) # max of first time deriv of pi
h5write(fileMax, "phidthdthmaxs",phidthdthmaxval[:]) # max of second angular deriv of phi
