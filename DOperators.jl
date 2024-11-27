########################
#Finite difference operators
########################
#Author: Alejandro Cardenas-Avendano

#First derivatives

#Fourth order Center Spatial Derivatives
const cd1y=[[-2,-1,1,2],[0,0,0,0]];
const cd1z=[[0,0,0,0],[-2,-1,1,2]];

const cd1yB1=[[0,-1,1,2],[0,0,0,0]];
const cd1yB2=[[2,1,1,2],[0,0,0,0]];

function cd_diff(vM,i2,i3,h,dc)
    return (vM[i2+dc[1][1],i3+dc[2][1]]-8.0*vM[i2+dc[1][2],i3+dc[2][2]]+8.0*vM[i2+dc[1][3],i3+dc[2][3]]-vM[i2+dc[1][4],i3+dc[2][4]])/(12.0*h)
end;

function cd_diffodd(vM,i2,i3,h,dc)
  return (-vM[i2+dc[1][1],i3+dc[2][1]]-8.0*vM[i2+dc[1][2],i3+dc[2][2]]+8.0*vM[i2+dc[1][3],i3+dc[2][3]]-vM[i2+dc[1][4],i3+dc[2][4]])/(12.0*h)
end;

#Sixth order Center Spatial Derivatives
#const cd1x6=[[-3,-2,-1,1,2,3],[0,0,0,0,0,0],[0,0,0,0,0,0]];
const cd1y6=[[-3,-2,-1,1,2,3],[0,0,0,0,0,0]];
const cd1z6=[[0,0,0,0,0,0],[-3,-2,-1,1,2,3]];

const cd1yB16=[[1,0,-1,1,2,3],[0,0,0,0,0,0]];
const cd1yB26=[[3,2,1,1,2,3],[0,0,0,0,0,0]];
const cd1yB36=[[-1,-2,-1,1,2,3],[0,0,0,0,0,0]];

function cd_diff6(vM,i2,i3,h,dc)
    return (-vM[i2+dc[1][1],i3+dc[2][1]]+9.0*vM[i2+dc[1][2],i3+dc[2][2]]-45.0*vM[i2+dc[1][3],i3+dc[2][3]]+45.0*vM[i2+dc[1][4],i3+dc[2][4]]-9.0*vM[i2+dc[1][5],i3+dc[2][5]]+vM[i2+dc[1][6],i3+dc[2][6]])/(60.0*h)
end;

#Fourth order Forward Spatial Derivatives
#const fw1x=[[0,1,2,3,4],[0,0,0,0,0],[0,0,0,0,0]];
const fw1y=[[0,1,2,3,4],[0,0,0,0,0]];
const fw1z=[[0,0,0,0,0],[0,1,2,3,4]];

function fw_diff(vM,i2,i3,h,dc)
  return (-25.0*vM[i2+dc[1][1],i3+dc[2][1]]+48.0*vM[i2+dc[1][2],i3+dc[2][2]]-36.0*vM[i2+dc[1][3],i3+dc[2][3]]+16.0*vM[i2+dc[1][4],i3+dc[2][4]]-3.0*vM[i2+dc[1][5],i3+dc[2][5]])/(12.0*h)
end;

#Fourth order Backward Spatial Derivatives
#const bw1x=[[-4,-3,-2,-1,0],[0,0,0,0,0],[0,0,0,0,0]];
const bw1y=[[-4,-3,-2,-1,0],[0,0,0,0,0]];
const bw1z=[[0,0,0,0,0],[-4,-3,-2,-1,0]];

function bw_diff(vM,i2,i3,h,dc)
  return (3.0*vM[i2+dc[1][1],i3+dc[2][1]]-16.0*vM[i2+dc[1][2],i3+dc[2][2]]+36.0*vM[i2+dc[1][3],i3+dc[2][3]]-48.0*vM[i2+dc[1][4],i3+dc[2][4]]+25.0*vM[i2+dc[1][5],i3+dc[2][5]])/(12.0*h)
end;

#Cartoon Neumann Boundary Conditions

const fw1yC4=[[1,2,3,4],[0,0,0,0]];

function fw_diffC4(vM,i2,i3,dc)
  return (48.0*vM[i2+dc[1][1],i3+dc[2][1]]-36.0*vM[i2+dc[1][2],i3+dc[2][2]]+16.0*vM[i2+dc[1][3],i3+dc[2][3]]-3.0*vM[i2+dc[1][4],i3+dc[2][4]])/(25.0)
end;

const fw1yC6=[[1,2,3,4,5,6],[0,0,0,0,0,0]];

function fw_diffC6(vM,i2,i3,dc)
  return (360.0*vM[i2+dc[1][1],i3+dc[2][1]]-450.0*vM[i2+dc[1][2],i3+dc[2][2]]+400.0*vM[i2+dc[1][3],i3+dc[2][3]]-225.0*vM[i2+dc[1][4],i3+dc[2][4]]+72.0*vM[i2+dc[1][5],i3+dc[2][5]]-10.0*vM[i2+dc[1][6],i3+dc[2][6]])/(147.0)
end;

const fw1yC8=[[1,2,3,4,5,6,7,8],[0,0,0,0,0,0,0,0]];

function fw_diffC8(vM,i2,i3,dc)
  return (6720.0*vM[i2+dc[1][1],i3+dc[2][1]]-11760.0*vM[i2+dc[1][2],i3+dc[2][2]]+15680.0*vM[i2+dc[1][3],i3+dc[2][3]]-14700.0*vM[i2+dc[1][4],i3+dc[2][4]]+9408.0*vM[i2+dc[1][5],i3+dc[2][5]]-3920.0*vM[i2+dc[1][6],i3+dc[2][6]]+960.0*vM[i2+dc[1][7],i3+dc[2][7]]-105.0*vM[i2+dc[1][8],i3+dc[2][8]])/(2283.0)
end;

####################
#Sixth derivatives
####################

#Kreiss-Oliger Dissipation 
#https://einsteintoolkit.org/thornguide/CactusNumerical/Dissipation/documentation.html

# #=
#Second order Center Spatial Derivatives
#const cd6x=[[-3,-2,-1,0,1,2,3],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]];
const cd6y=[[-3,-2,-1,0,1,2,3],[0,0,0,0,0,0,0]];
const cd6z=[[0,0,0,0,0,0,0],[-3,-2,-1,0,1,2,3]];

const cd6yL1=[[-1,-2,-1,0,1,2,3],[0,0,0,0,0,0,0]];
const cd6yL2=[[1,0,-1,0,1,2,3],[0,0,0,0,0,0,0]];
const cd6yL3=[[3,2,1,0,1,2,3],[0,0,0,0,0,0,0]];

const cd6zL1=[[0,0,0,0,0,0,0],[-1,-2,-1,0,1,2,3]];
const cd6zL2=[[0,0,0,0,0,0,0],[1,0,-1,0,1,2,3]];
#const cd6zL3=[[0,0,0,0,0,0,0],[3,2,1,0,1,2,3]];

const cd6yR1=[[-3,-2,-1,0,1,2,1],[0,0,0,0,0,0,0]];
const cd6yR2=[[-3,-2,-1,0,1,0,-1],[0,0,0,0,0,0,0]];

const cd6zR1=[[0,0,0,0,0,0,0],[-3,-2,-1,0,1,2,1]];
const cd6zR2=[[0,0,0,0,0,0,0],[-3,-2,-1,0,1,0,-1]];
#const cd6zR3=[[0,0,0,0,0,0,0],[-3,-2,-1,0,1,2,3]];

function cd_diss6(vM,i2,i3,dc)
  return epsKO/(64.0)*(vM[i2+dc[1][1],i3+dc[2][1]] - 6.0*vM[i2+dc[1][2],i3+dc[2][2]]+15.0*vM[i2+dc[1][3],i3+dc[2][3]]-20.0*vM[i2+dc[1][4],i3+dc[2][4]]+15.0*vM[i2+dc[1][5],i3+dc[2][5]]-6.0*vM[i2+dc[1][6],i3+dc[2][6]]+vM[i2+dc[1][7],i3+dc[2][7]])
end;
