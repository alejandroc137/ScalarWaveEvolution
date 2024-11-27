########################
#Initial Conditions
########################
#Author: Alejandro Cardenas-Avendano

function arg(y,z)
    if z==0.0
        0.0
    else
        tan(pi*z/2.0)/sqrt((tan(pi*y/2.0))^2.0+(tan(pi*z/2.0))^2.0)
    end
end

function Phival(z,y)
    phivalue=0
    if r0<=sqrt(tan((pi*y)/2.0)^2.0 + tan((pi*z)/2.0)^2.0)<=r1
        for ell=1:size(ells)[1]
            phivalue+=GSL.sf_legendre_sphPlm(ells[ell],0,arg(y,z))*(Am*((-r0 + sqrt(tan((pi*y)/2.0)^2.0 + tan((pi*z)/2.0)^2.0))^polylogexp*(-r1 + sqrt(tan((pi*y)/2.0)^2.0 + tan((pi*z)/2.0)^2.0))^polylogexp)/(r0-r1)^(2.0*polylogexp))
        end
    end
    return phivalue
end

function ICs!(phi_M1n,pi_M1n,gttn,gxxn,gxyn,gxzn,gyyn,gyzn,gzzn,sqrtmingn)
    @inbounds Threads.@threads for iy=1:Ny0
        @inbounds for iz=1:Nz0
                
                gttn[iy,iz]=-((a + b*(tan((pi*0.0)/2.0)^2.0 + tan((pi*ys[iy])/2.0)^2.0 + tan((pi*zs[iz])/2.0)^2.0)^2.0)/(a + (tan((pi*0.0)/2.0)^2.0 + tan((pi*ys[iy])/2.0)^2.0 + tan((pi*zs[iz])/2.0)^2.0)*(-1.0 + b*(tan((pi*0.0)/2.0)^2.0 + tan((pi*ys[iy])/2.0)^2.0 + tan((pi*zs[iz])/2.0)^2.0))));

                gxxn[iy,iz]=(4.0*cos((pi*0.0)/2.0)^4.0*(a + b*tan((pi*0.0)/2.0)^4.0 + b*tan((pi*ys[iy])/2.0)^4.0 + 2.0*b*tan((pi*ys[iy])/2.0)^2.0*tan((pi*zs[iz])/2.0)^2.0 + b*tan((pi*zs[iz])/2.0)^4.0 + tan((pi*0.0)/2.0)^2.0*(-1.0 + 2.0*b*tan((pi*ys[iy])/2.0)^2.0 + 2.0*b*tan((pi*zs[iz])/2.0)^2.0)))/(pi^2.0*(a + b*tan((pi*0.0)/2.0)^4.0 + b*tan((pi*ys[iy])/2.0)^4.0 + 2.0*b*tan((pi*ys[iy])/2.0)^2.0*tan((pi*zs[iz])/2.0)^2.0 + b*tan((pi*zs[iz])/2.0)^4.0 + 2.0*b*tan((pi*0.0)/2.0)^2.0*(tan((pi*ys[iy])/2.0)^2.0 + tan((pi*zs[iz])/2.0)^2.0)));

                gxyn[iy,iz]=-((sin(pi*0.0)*sin(pi*ys[iy]))/(pi^2.0*(a + b*tan((pi*0.0)/2.0)^4.0 + b*tan((pi*ys[iy])/2.0)^4.0 + 2.0*b*tan((pi*ys[iy])/2.0)^2.0*tan((pi*zs[iz])/2.0)^2.0 + b*tan((pi*zs[iz])/2.0)^4.0 + 2.0*b*tan((pi*0.0)/2.0)^2.0*(tan((pi*ys[iy])/2.0)^2.0 + tan((pi*zs[iz])/2.0)^2.0))));

                gxzn[iy,iz]=-((sin(pi*0.0)*sin(pi*zs[iz]))/(pi^2.0*(a + b*tan((pi*0.0)/2.0)^4.0 + b*tan((pi*ys[iy])/2.0)^4.0 + 2.0*b*tan((pi*ys[iy])/2.0)^2.0*tan((pi*zs[iz])/2.0)^2.0 + b*tan((pi*zs[iz])/2.0)^4.0 + 2.0*b*tan((pi*0.0)/2.0)^2.0*(tan((pi*ys[iy])/2.0)^2.0 + tan((pi*zs[iz])/2.0)^2.0))));

                gyyn[iy,iz]=(4.0*cos((pi*ys[iy])/2.0)^4.0*(a + b*tan((pi*0.0)/2.0)^4.0 + b*tan((pi*ys[iy])/2.0)^4.0 + b*tan((pi*zs[iz])/2.0)^4.0 + 2.0*b*tan((pi*0.0)/2.0)^2.0*(tan((pi*ys[iy])/2.0)^2.0 + tan((pi*zs[iz])/2.0)^2.0) + tan((pi*ys[iy])/2.0)^2.0*(-1.0 + 2.0*b*tan((pi*zs[iz])/2.0)^2.0)))/(pi^2.0*(a + b*tan((pi*0.0)/2.0)^4.0 + b*tan((pi*ys[iy])/2.0)^4.0 + 2.0*b*tan((pi*ys[iy])/2.0)^2.0*tan((pi*zs[iz])/2.0)^2.0 + b*tan((pi*zs[iz])/2.0)^4.0 + 2.0*b*tan((pi*0.0)/2.0)^2.0*(tan((pi*ys[iy])/2.0)^2.0 + tan((pi*zs[iz])/2.0)^2.0)));

                gyzn[iy,iz]=-((sin(pi*ys[iy])*sin(pi*zs[iz]))/(pi^2.0*(a + b*tan((pi*0.0)/2.0)^4.0 + b*tan((pi*ys[iy])/2.0)^4.0 + 2.0*b*tan((pi*ys[iy])/2.0)^2.0*tan((pi*zs[iz])/2.0)^2.0 + b*tan((pi*zs[iz])/2.0)^4.0 + 2.0*b*tan((pi*0.0)/2.0)^2.0*(tan((pi*ys[iy])/2.0)^2.0 + tan((pi*zs[iz])/2.0)^2.0))));

                gzzn[iy,iz]=(4.0*cos((pi*zs[iz])/2.0)^4.0*(a + b*tan((pi*0.0)/2.0)^4.0 + b*tan((pi*ys[iy])/2.0)^4.0 - tan((pi*zs[iz])/2.0)^2.0 + 2.0*b*tan((pi*ys[iy])/2.0)^2.0*tan((pi*zs[iz])/2.0)^2.0 + b*tan((pi*zs[iz])/2.0)^4.0 + 2.0*b*tan((pi*0.0)/2.0)^2.0*(tan((pi*ys[iy])/2.0)^2.0 + tan((pi*zs[iz])/2.0)^2.0)))/(pi^2.0*(a + b*tan((pi*0.0)/2.0)^4.0 + b*tan((pi*ys[iy])/2.0)^4.0 + 2.0*b*tan((pi*ys[iy])/2.0)^2.0*tan((pi*zs[iz])/2.0)^2.0 + b*tan((pi*zs[iz])/2.0)^4.0 + 2.0*b*tan((pi*0.0)/2.0)^2.0*(tan((pi*ys[iy])/2.0)^2.0 + tan((pi*zs[iz])/2.0)^2.0)));
                            
                sqrtmingn[iy,iz]=(8.0*(cos((pi*0.0)/2.0)^2.0*cos((pi*ys[iy])/2.0)^2.0*cos((pi*zs[iz])/2.0)^2.0)/pi^3.0)^(-1.0);
        end
    end
    if ((Ti==0.0) & (ICtype!="File")) 
        println("======================")
        println("Using new IC with Polylog conditions")
        println("Exponent $(polylogexp) and ells $(ells)")
        println("======================")
        @inbounds Threads.@threads for iy=1:Ny0
            @inbounds for iz=1:Nz0
                    #Polylog with exponent given by polylogexp and ells for the spherical
                    phi_M1n[iy,iz]=Phival(zs[iz],ys[iy])
            end
        end

    elseif (Ti!=0.0) 
        println("Using ICs from file snapshot at time $(Ti)")
        fid = h5open(location*"4DWave_$(filename)_$(Am)_$(Ti)_$(dy)_$(lambda)_$(epsKO).h5");
        phi_M1n[:,:]=read(fid["phi"]);
        pi_M1n[:,:]=read(fid["pi"]);
    end
end;