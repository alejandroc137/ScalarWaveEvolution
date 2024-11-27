#This code solves the non-linear wave equation in 2+1. 
#It implements a RK4 to evolve in time the system, 
#and fourth order finite differences for the spatial derivatives
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
            default = 4000.0
    end
    return parse_args(s)
end

parsed_args = parse_commandline()

spacetime="NLMWP"
ICtype="Analytical"

orderFD=convert(Int64,4)


#For the piecewise poly
const Am=parsed_args["amp"];
const r0=0.02;
const r1=0.6;
const ells=[1,2]
const polylogexp=4.0

# Resolution of the grids
const dy=parsed_args["res"];
const dz=parsed_args["res"];

# Courant factor: dt/dx
const lambda = 0.5;
# Kreiss--Oliger dissipation
const epsKO=0.3;
# Initial Time
const Ti= parsed_args["ti"];
# Final Time
const Tf= parsed_args["tf"];

println("Initial time  ", Ti)
println("Final time  ", Tf)

ymin       =  0.0;
ymax       =  1.0;

zmin       =  -1.0;
zmax       =  1.0;

# Define the folder name
folder_name = "Results"

# Check if the folder exists
if !isdir(folder_name)
    # Create the folder if it doesn't exist
    mkdir(folder_name)
    println("Folder 'Results' created.")
else
    println("Folder 'Results' already exists.")
end

location="Results/" 

filename=spacetime

#For the metric potential
if spacetime=="Minkowski"
    const a=1e8;
    const b=1e8;
else
    const a=0.026;
    const b=11.20;
end

const dt=lambda*dy;

const Ny0= convert(Int64,(ymax-ymin)/dy+1);
const Nz0= convert(Int64,(zmax-zmin)/dz+1);


println("Allocating memory for the matrices of size $Ny0 x $Nz0");

const Nt0= convert(Int64,round((Tf-Ti)/dt));

const ts=collect(range(Ti, Tf, length=Nt0+1));

const roundfact=1e8

#Cadence to auxiliar quantites 
const dtRaux=0.05;
#Cadence to store the field and its time derivative
const dtRaux2=0.1;

const NtR=convert(Int64,round((Tf-Ti)/dtRaux));
#Cadence to save the max vals
const dtR=round(dtRaux*roundfact);
#Cadence to save images of the simulation
const dtRFS=round(dtRaux2*roundfact);

const ts_aux=round.(roundfact*ts);

const ys=collect(range(ymin, ymax, length=Ny0));
const zs=collect(range(zmin, zmax, length=Nz0));

const pidotmaxval=collect(range(Ti, Tf, length=Nt0+1));

const phidxdxmaxval=collect(range(Ti, Tf, length=Nt0+1));
const phidydymaxval=collect(range(Ti, Tf, length=Nt0+1));
const phidzdzmaxval=collect(range(Ti, Tf, length=Nt0+1));
const phidydzmaxval=collect(range(Ti, Tf, length=Nt0+1));

#Spherical Coordinates

const phidthdthmaxval=collect(range(Ti, Tf, length=Nt0+1));
const phidphidphimaxval=collect(range(Ti, Tf, length=Nt0+1));

boundenerg=parsed_args["energb"];

yvaltest=findall(x -> x == boundenerg, ys)[1];

zvaltest1=findall(x -> x == -boundenerg, zs)[1];
zvaltest2=findall(x -> x == boundenerg, zs)[1];

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

#Contravariant
const gtt=zeros(Float64, (Ny0,Nz0));

const gxx=zeros(Float64, (Ny0,Nz0));
const gxy=zeros(Float64, (Ny0,Nz0));
const gxz=zeros(Float64, (Ny0,Nz0));

const gyy=zeros(Float64, (Ny0,Nz0));
const gyz=zeros(Float64, (Ny0,Nz0));

const gzz=zeros(Float64, (Ny0,Nz0));

#Derivatives

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

fileInfo=location*"Info_Wave_$(filename)_$(Am)_$(dy)_$(lambda)_$(epsKO)_$(Ti)_$(Tf).txt"
fileMax=location*"Max_Wave_$(filename)_$(Am)_$(dy)_$(lambda)_$(epsKO)_$(Ti)_$(Tf).h5"

if isfile(fileInfo)
    rm(fileInfo)
    println("Info File Overwritten")
end

write(fileInfo," Simulation params \n\n Am=$(Am) \n r0=$(r0)\n r1=$(r1)\n ells=$(ells)\n polylogexp=$(polylogexp)\n a=$(a)\n b=$(b) " )

println("Evolving the wave equation");
@time simulation!(size(ts)[1],dtR,phi_M1,pi_M1,phi_M2,pi_M2,phidx_M,pidx_M,phidy_M,pidy_M,phidz_M,pidz_M,phidot_M,pidot_M,phidxdx_M,phidxdy_M,phidxdz_M,phidydy_M,phidydz_M,phidzdz_M)

if isfile(fileMax)
    rm(fileMax)
    println("File Maxs Overwritten")
end
h5write(fileMax, "ts",ts[:])
h5write(fileMax, "pidotmaxs",pidotmaxval[:])
h5write(fileMax, "phidthdthmaxs",phidthdthmaxval[:])
