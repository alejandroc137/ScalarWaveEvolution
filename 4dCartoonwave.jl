########################
#Discretization of the equations 
########################
#Author: Alejandro Cardenas-Avendano

#Metric derivatives 
function metricderivatives!(gtt_dxn,gxx_dxn,gxy_dxn,gxz_dxn,gyy_dxn,gyz_dxn,gzz_dxn,sqrtming_dxn,gtt_dyn,gxx_dyn,gxy_dyn,gxz_dyn,gyy_dyn,gyz_dyn,gzz_dyn,sqrtming_dyn,gtt_dzn,gxx_dzn,gxy_dzn,gxz_dzn,gyy_dzn,gyz_dzn,gzz_dzn,sqrtming_dzn)
    
    #X derivatives
    @inbounds Threads.@threads for iy=2:(Ny0-1)
        @inbounds for iz=2:(Nz0-1)
            gxx_dxn[iy,iz]=2.0*pi*csc(pi*ys[iy])*gxy[iy,iz]
            gxy_dxn[iy,iz]=(pi*cot(pi*ys[iy]/2.0)*(-cos(pi*ys[iy]/2.0)^4.0*gxx[iy,iz]+gyy[iy,iz]))/(1.0+cos(pi*ys[iy]))
            gxz_dxn[iy,iz]=pi*cot(pi*ys[iy]/2.0)*gyz[iy,iz]/(1.0+cos(pi*ys[iy]))
            gyy_dxn[iy,iz]=-pi*cos(pi*ys[iy]/2.0)^2.0*cot(pi*ys[iy]/2.0)*gxy[iy,iz]
            gyz_dxn[iy,iz]=-pi*cos(pi*ys[iy]/2.0)^2.0*cot(pi*ys[iy]/2.0)*gxz[iy,iz]/2.0
        end
    end

    #Y derivatives
    @inbounds Threads.@threads for iz=2:(Nz0-1)

        gtt_dyn[1,iz]=fw_diff(gtt,1,iz,dy,fw1y)
        gxx_dyn[1,iz]=fw_diff(gxx,1,iz,dy,fw1y)
        gxy_dyn[1,iz]=fw_diff(gxy,1,iz,dy,fw1y)
        gxz_dyn[1,iz]=fw_diff(gxz,1,iz,dy,fw1y)
        gyy_dyn[1,iz]=fw_diff(gyy,1,iz,dy,fw1y)
        gyz_dyn[1,iz]=fw_diff(gyz,1,iz,dy,fw1y)
        gzz_dyn[1,iz]=fw_diff(gzz,1,iz,dy,fw1y)
        sqrtming_dyn[1,iz]=fw_diff(sqrtming,1,iz,dy,fw1y)

        gtt_dyn[2,iz]=fw_diff(gtt,2,iz,dy,fw1y)
        gxx_dyn[2,iz]=fw_diff(gxx,2,iz,dy,fw1y)
        gxy_dyn[2,iz]=fw_diff(gxy,2,iz,dy,fw1y)
        gxz_dyn[2,iz]=fw_diff(gxz,2,iz,dy,fw1y)
        gyy_dyn[2,iz]=fw_diff(gyy,2,iz,dy,fw1y)
        gyz_dyn[2,iz]=fw_diff(gyz,2,iz,dy,fw1y)
        gzz_dyn[2,iz]=fw_diff(gzz,2,iz,dy,fw1y)
        sqrtming_dyn[2,iz]=fw_diff(sqrtming,2,iz,dy,fw1y)

        gtt_dyn[Ny0-1,iz]=bw_diff(gtt,Ny0-1,iz,dy,bw1y)
        gxx_dyn[Ny0-1,iz]=bw_diff(gxx,Ny0-1,iz,dy,bw1y)
        gxy_dyn[Ny0-1,iz]=bw_diff(gxy,Ny0-1,iz,dy,bw1y)
        gxz_dyn[Ny0-1,iz]=bw_diff(gxz,Ny0-1,iz,dy,bw1y)
        gyy_dyn[Ny0-1,iz]=bw_diff(gyy,Ny0-1,iz,dy,bw1y)
        gyz_dyn[Ny0-1,iz]=bw_diff(gyz,Ny0-1,iz,dy,bw1y)
        gzz_dyn[Ny0-1,iz]=bw_diff(gzz,Ny0-1,iz,dy,bw1y)
        sqrtming_dyn[Ny0-1,iz]=bw_diff(sqrtming,Ny0-1,iz,dy,bw1y)

        gtt_dyn[Ny0-2,iz]=bw_diff(gtt,Ny0-2,iz,dy,bw1y)
        gxx_dyn[Ny0-2,iz]=bw_diff(gxx,Ny0-2,iz,dy,bw1y)
        gxy_dyn[Ny0-2,iz]=bw_diff(gxy,Ny0-2,iz,dy,bw1y)
        gxz_dyn[Ny0-2,iz]=bw_diff(gxz,Ny0-2,iz,dy,bw1y)
        gyy_dyn[Ny0-2,iz]=bw_diff(gyy,Ny0-2,iz,dy,bw1y)
        gyz_dyn[Ny0-2,iz]=bw_diff(gyz,Ny0-2,iz,dy,bw1y)
        gzz_dyn[Ny0-2,iz]=bw_diff(gzz,Ny0-2,iz,dy,bw1y)
        sqrtming_dyn[Ny0-2,iz]=bw_diff(sqrtming,Ny0-2,iz,dy,bw1y)

    end

    @inbounds Threads.@threads for iy=3:(Ny0-3)
        @inbounds for iz=2:(Nz0-1)
            gtt_dyn[iy,iz]=cd_diff(gtt,iy,iz,dy,cd1y)
            gxx_dyn[iy,iz]=cd_diff(gxx,iy,iz,dy,cd1y)
            gxy_dyn[iy,iz]=cd_diff(gxy,iy,iz,dy,cd1y)
            gxz_dyn[iy,iz]=cd_diff(gxz,iy,iz,dy,cd1y)
            gyy_dyn[iy,iz]=cd_diff(gyy,iy,iz,dy,cd1y)
            gyz_dyn[iy,iz]=cd_diff(gyz,iy,iz,dy,cd1y)
            gzz_dyn[iy,iz]=cd_diff(gzz,iy,iz,dy,cd1y)
            sqrtming_dyn[iy,iz]=cd_diff(sqrtming,iy,iz,dy,cd1y)
        end
    end

    #Z derivatives
    @inbounds Threads.@threads for iy=1:(Ny0-1)

        gtt_dzn[iy,2]=fw_diff(gtt,iy,2,dz,fw1z)
        gxx_dzn[iy,2]=fw_diff(gxx,iy,2,dz,fw1z)
        gxy_dzn[iy,2]=fw_diff(gxy,iy,2,dz,fw1z)
        gxz_dzn[iy,2]=fw_diff(gxz,iy,2,dz,fw1z)
        gyy_dzn[iy,2]=fw_diff(gyy,iy,2,dz,fw1z)
        gyz_dzn[iy,2]=fw_diff(gyz,iy,2,dz,fw1z)
        gzz_dzn[iy,2]=fw_diff(gzz,iy,2,dz,fw1z)
        sqrtming_dzn[iy,2]=fw_diff(sqrtming,iy,2,dz,fw1z)

        gtt_dzn[iy,3]=fw_diff(gtt,iy,3,dz,fw1z)
        gxx_dzn[iy,3]=fw_diff(gxx,iy,3,dz,fw1z)
        gxy_dzn[iy,3]=fw_diff(gxy,iy,3,dz,fw1z)
        gxz_dzn[iy,3]=fw_diff(gxz,iy,3,dz,fw1z)
        gyy_dzn[iy,3]=fw_diff(gyy,iy,3,dz,fw1z)
        gyz_dzn[iy,3]=fw_diff(gyz,iy,3,dz,fw1z)
        gzz_dzn[iy,3]=fw_diff(gzz,iy,3,dz,fw1z)
        sqrtming_dzn[iy,3]=fw_diff(sqrtming,iy,3,dz,fw1z)

        gtt_dzn[iy,Nz0-1]=bw_diff(gtt,iy,Nz0-1,dz,bw1z)
        gxx_dzn[iy,Nz0-1]=bw_diff(gxx,iy,Nz0-1,dz,bw1z)
        gxy_dzn[iy,Nz0-1]=bw_diff(gxy,iy,Nz0-1,dz,bw1z)
        gxz_dzn[iy,Nz0-1]=bw_diff(gxz,iy,Nz0-1,dz,bw1z)
        gyy_dzn[iy,Nz0-1]=bw_diff(gyy,iy,Nz0-1,dz,bw1z)
        gyz_dzn[iy,Nz0-1]=bw_diff(gyz,iy,Nz0-1,dz,bw1z)
        gzz_dzn[iy,Nz0-1]=bw_diff(gzz,iy,Nz0-1,dz,bw1z)
        sqrtming_dzn[iy,Nz0-1]=bw_diff(sqrtming,iy,Nz0-1,dz,bw1z)

        gtt_dzn[iy,Nz0-2]=bw_diff(gtt,iy,Nz0-2,dz,bw1z)
        gxx_dzn[iy,Nz0-2]=bw_diff(gxx,iy,Nz0-2,dz,bw1z)
        gxy_dzn[iy,Nz0-2]=bw_diff(gxy,iy,Nz0-2,dz,bw1z)
        gxz_dzn[iy,Nz0-2]=bw_diff(gxz,iy,Nz0-2,dz,bw1z)
        gyy_dzn[iy,Nz0-2]=bw_diff(gyy,iy,Nz0-2,dz,bw1z)
        gyz_dzn[iy,Nz0-2]=bw_diff(gyz,iy,Nz0-2,dz,bw1z)
        gzz_dzn[iy,Nz0-2]=bw_diff(gzz,iy,Nz0-2,dz,bw1z)
        sqrtming_dzn[iy,Nz0-2]=bw_diff(sqrtming,iy,Nz0-2,dz,bw1z)
    end

    @inbounds Threads.@threads for iy=1:(Ny0-1)
        @inbounds for iz=4:(Nz0-3)
            gtt_dzn[iy,iz]=cd_diff(gtt,iy,iz,dz,cd1z)
            gxx_dzn[iy,iz]=cd_diff(gxx,iy,iz,dz,cd1z)
            gxy_dzn[iy,iz]=cd_diff(gxy,iy,iz,dz,cd1z)
            gxz_dzn[iy,iz]=cd_diff(gxz,iy,iz,dz,cd1z)
            gyy_dzn[iy,iz]=cd_diff(gyy,iy,iz,dz,cd1z)
            gyz_dzn[iy,iz]=cd_diff(gyz,iy,iz,dz,cd1z)
            gzz_dzn[iy,iz]=cd_diff(gzz,iy,iz,dz,cd1z)
            sqrtming_dzn[iy,iz]=cd_diff(sqrtming,iy,iz,dz,cd1z)
        end
    end  
end

function RHS_Eq!(phi_Mn,pi_Mn,phidx_Mn,pidx_Mn,phidy_Mn,pidy_Mn,phidz_Mn,pidz_Mn,phidot_Mn,pidot_Mn,phidxdx_Mn,phidxdy_Mn,phidxdz_Mn,phidydy_Mn,phidydz_Mn,phidzdz_Mn)

    #First derivatives
    #Y Derivatives
    @inbounds Threads.@threads for iz=2:(Nz0-1)

        #phidy_Mn[1,iz]=cd_diff(phi_Mn,1,iz,dy,cd1yB2)
        phidy_Mn[2,iz]=cd_diff(phi_Mn,2,iz,dy,cd1yB1)

        phidy_Mn[Ny0-1,iz]=bw_diff(phi_Mn,Ny0-1,iz,dy,bw1y)
        phidy_Mn[Ny0-2,iz]=bw_diff(phi_Mn,Ny0-2,iz,dy,bw1y)

    end

    @inbounds Threads.@threads for iy=3:(Ny0-3)
        @inbounds for iz=2:(Nz0-1)
            phidy_Mn[iy,iz]=cd_diff(phi_Mn,iy,iz,dy,cd1y)
            #pidy_Mn[iy,iz]=cd_diff(pi_Mn,1,iy,iz,dy,cd1y)
        end
    end

    #Z Derivatives
    @inbounds Threads.@threads for iy=1:(Ny0-1)

        phidz_Mn[iy,2]=fw_diff(phi_Mn,iy,2,dz,fw1z)
        #pidz_Mn[iy,2]=fw_diff(pi_Mn,iy,2,dz,fw1z)

        phidz_Mn[iy,3]=fw_diff(phi_Mn,iy,3,dz,fw1z)
        #pidz_Mn[iy,3]=fw_diff(pi_Mn,iy,3,dz,fw1z)

        phidz_Mn[iy,Nz0-1]=bw_diff(phi_Mn,iy,Nz0-1,dz,bw1z)
        #pidz_Mn[iy,Nz0-1]=bw_diff(pi_Mn,iy,Nz0-1,dz,bw1z)

        phidz_Mn[iy,Nz0-2]=bw_diff(phi_Mn,iy,Nz0-2,dz,bw1z)
        #pidz_Mn[iy,Nz0-2]=bw_diff(pi_Mn,iy,Nz0-2,dz,bw1z)

    end

    @inbounds Threads.@threads for iy=1:(Ny0-1)
        @inbounds for iz=4:(Nz0-3)
            phidz_Mn[iy,iz]=cd_diff(phi_Mn,iy,iz,dz,cd1z)
            #pidz_Mn[iy,iz]=cd_diff(pi_Mn,1,iy,iz,dz,cd1z)
        end
    end  

    ######################
    ##Second derivatives##
    ######################

    #Y Derivatives

    @inbounds Threads.@threads for iz=2:(Nz0-1)

        phidxdx_Mn[2,iz]=(pi/16.0)*(csc(pi*ys[2]/2.0)^4.0*sin(pi*ys[2])^3.0)*phidy_Mn[2,iz]        
        phidxdx_Mn[Ny0-1,iz]=(pi/16.0)*(csc(pi*ys[Ny0-1]/2.0)^4.0*sin(pi*ys[Ny0-1])^3.0)*phidy_Mn[Ny0-1,iz]
        phidxdx_Mn[Ny0-2,iz]=(pi/16.0)*(csc(pi*ys[Ny0-2]/2.0)^4.0*sin(pi*ys[Ny0-2])^3.0)*phidy_Mn[Ny0-2,iz]
                             
        phidydy_Mn[2,iz]=fw_diff(phidy_Mn,2,iz,dy,fw1y)
        phidydy_Mn[Ny0-1,iz]=bw_diff(phidy_Mn,Ny0-1,iz,dy,bw1y)
        phidydy_Mn[Ny0-2,iz]=bw_diff(phidy_Mn,Ny0-2,iz,dy,bw1y)

        phidydz_Mn[2,iz]=fw_diff(phidz_Mn,2,iz,dy,fw1y)
        phidydz_Mn[Ny0-1,iz]=bw_diff(phidz_Mn,Ny0-1,iz,dy,bw1y)
        phidydz_Mn[Ny0-2,iz]=bw_diff(phidz_Mn,Ny0-2,iz,dy,bw1y)
    end

    @inbounds Threads.@threads for iy=3:(Ny0-3)
        @inbounds for iz=2:Nz0-1
            phidxdx_Mn[iy,iz]=(pi/16.0)*(csc(pi*ys[iy]/2.0)^4.0*sin(pi*ys[iy])^3.0)*phidy_Mn[iy,iz]
            phidydy_Mn[iy,iz]=cd_diff(phidy_Mn,iy,iz,dy,cd1y)
            phidydz_Mn[iy,iz]=cd_diff(phidz_Mn,iy,iz,dy,cd1y)
        end
    end

    #Z Derivatives
    @inbounds Threads.@threads for iy=2:(Ny0-1)
        @inbounds for iz=4:(Nz0-3)
            phidzdz_Mn[iy,iz]=cd_diff(phidz_Mn,iy,iz,dz,cd1z)
        end
    end

    @inbounds Threads.@threads for iy=2:(Ny0-1)
        phidzdz_Mn[iy,2]=fw_diff(phidz_Mn,iy,2,dz,fw1z)
        phidzdz_Mn[iy,3]=fw_diff(phidz_Mn,iy,3,dz,fw1z)

        phidzdz_Mn[iy,Nz0-1]=bw_diff(phidz_Mn,iy,Nz0-1,dz,bw1z)
        phidzdz_Mn[iy,Nz0-2]=bw_diff(phidz_Mn,iy,Nz0-2,dz,bw1z)
    end

    #Evolution equations
    @inbounds Threads.@threads for iy=2:(Ny0-1)
        @inbounds for iz=2:(Nz0-1)

            phidot_Mn[iy,iz]=-pi_Mn[iy,iz]/gtt[iy,iz];

            #Linear Wave Equation
            #pidot_Mn[iy,iz]=(                                                 gyz[iy,iz]*sqrtming_dy[iy,iz]*phidz_Mn[iy,iz] + gxz[iy,iz]*sqrtming_dx[iy,iz]*phidz_Mn[iy,iz] + gzz[iy,iz]*(sqrtming_dz[iy,iz]*phidz_Mn[iy,iz] + sqrtming[iy,iz]*phidzdz_Mn[iy,iz])+ gyz[iy,iz]*sqrtming_dz[iy,iz]*phidy_Mn[iy,iz] + gyy[iy,iz]*sqrtming_dy[iy,iz]*phidy_Mn[iy,iz] + gxy[iy,iz]*sqrtming_dx[iy,iz]*phidy_Mn[iy,iz] + (gxz[iy,iz]*sqrtming_dz[iy,iz] + gxy[iy,iz]*sqrtming_dy[iy,iz] + gxx[iy,iz]*sqrtming_dx[iy,iz])*phidx_Mn[iy,iz] + sqrtming[iy,iz]*((gzz_dz[iy,iz] + gyz_dy[iy,iz] + gxz_dx[iy,iz])*phidz_Mn[iy,iz] + (gyz_dz[iy,iz] + gyy_dy[iy,iz] + gxy_dx[iy,iz])*phidy_Mn[iy,iz] + 2.0*gyz[iy,iz]*phidydz_Mn[iy,iz] + gyy[iy,iz]*phidydy_Mn[iy,iz] + (gxz_dz[iy,iz] + gxy_dy[iy,iz] + gxx_dx[iy,iz])*phidx_Mn[iy,iz]  + 2.0*gxz[iy,iz]*phidxdz_Mn[iy,iz] + 2.0*gxy[iy,iz]*phidxdy_Mn[iy,iz] + gxx[iy,iz]*phidxdx_Mn[iy,iz]))/sqrtming[iy,iz]
            
            #Non-Linear Wave Equation 
            pidot_Mn[iy,iz]=(-(sqrtming[iy,iz]*phi_Mn[iy,iz]^3.0)    + gyz[iy,iz]*sqrtming_dy[iy,iz]*phidz_Mn[iy,iz] + gxz[iy,iz]*sqrtming_dx[iy,iz]*phidz_Mn[iy,iz] + gzz[iy,iz]*(sqrtming_dz[iy,iz]*phidz_Mn[iy,iz] + sqrtming[iy,iz]*phidzdz_Mn[iy,iz])+ gyz[iy,iz]*sqrtming_dz[iy,iz]*phidy_Mn[iy,iz] + gyy[iy,iz]*sqrtming_dy[iy,iz]*phidy_Mn[iy,iz] + gxy[iy,iz]*sqrtming_dx[iy,iz]*phidy_Mn[iy,iz] + (gxz[iy,iz]*sqrtming_dz[iy,iz] + gxy[iy,iz]*sqrtming_dy[iy,iz] + gxx[iy,iz]*sqrtming_dx[iy,iz])*phidx_Mn[iy,iz] + sqrtming[iy,iz]*((gzz_dz[iy,iz] + gyz_dy[iy,iz] + gxz_dx[iy,iz])*phidz_Mn[iy,iz] + (gyz_dz[iy,iz] + gyy_dy[iy,iz] + gxy_dx[iy,iz])*phidy_Mn[iy,iz] + 2.0*gyz[iy,iz]*phidydz_Mn[iy,iz] + gyy[iy,iz]*phidydy_Mn[iy,iz] + (gxz_dz[iy,iz] + gxy_dy[iy,iz] + gxx_dx[iy,iz])*phidx_Mn[iy,iz]  + 2.0*gxz[iy,iz]*phidxdz_Mn[iy,iz] + 2.0*gxy[iy,iz]*phidxdy_Mn[iy,iz] + gxx[iy,iz]*phidxdx_Mn[iy,iz]))/sqrtming[iy,iz]
            
        end
    end
end

function spatial!(f,dyf,dzf)

    @inbounds Threads.@threads for iz=2:(Nz0-1)

        dyf[2,iz]=cd_diff(f,2,iz,dy,cd1yB1)

        dyf[Ny0-1,iz]=bw_diff(f,Ny0-1,iz,dy,bw1y)
        dyf[Ny0-2,iz]=bw_diff(f,Ny0-2,iz,dy,bw1y)

    end

    @inbounds Threads.@threads for iy=3:(Ny0-3)
        @inbounds for iz=2:(Nz0-1)
            dyf[iy,iz]=cd_diff(f,iy,iz,dy,cd1y)
        end
    end

    #Z Derivatives
    @inbounds Threads.@threads for iy=1:(Ny0-1)

        dzf[iy,2]=fw_diff(f,iy,2,dz,fw1z)

        dzf[iy,3]=fw_diff(f,iy,3,dz,fw1z)

        dzf[iy,Nz0-1]=bw_diff(f,iy,Nz0-1,dz,bw1z)

        dzf[iy,Nz0-2]=bw_diff(f,iy,Nz0-2,dz,bw1z)

    end

    @inbounds Threads.@threads for iy=1:(Ny0-1)
        @inbounds for iz=4:(Nz0-3)
            dzf[iy,iz]=cd_diff(f,iy,iz,dz,cd1z)
        end
    end

end

function energyC(Qdx,Qdy,Qdz,Qdot)
    nrg=0.0;
    @inbounds for iy=2:(Ny0-1)
        @inbounds for iz=2:(Nz0-1)
                #Linear

                fact=dy*dz*1.0/sqrt(1.0 - (tan((pi*ys[iy])/2.0)^2.0 + tan((pi*zs[iz])/2.0)^2.0)/(a + b*(tan((pi*ys[iy])/2.0)^2.0 + tan((pi*zs[iz])/2.0)^2.0)^2.0))*(2.0*pi*tan(pi/2.0*ys[iy]))*(pi^2.0/4.0*sec(pi*ys[iy]/2.0)^2.0*sec(pi*zs[iz]/2.0)^2.0)
                
                nrg+=fact*((gzz[iy,iz]*Qdz[iy,iz]^2.0 + 2.0*gyz[iy,iz]*Qdz[iy,iz]*Qdy[iy,iz] + gyy[iy,iz]*Qdy[iy,iz]^2.0 + 2.0*(gxz[iy,iz]*phidz_M[iy,iz] + gxy[iy,iz]*Qdy[iy,iz])*Qdx[iy,iz] + gxx[iy,iz]*Qdx[iy,iz]^2.0 - gtt[iy,iz]*Qdot[iy,iz]^2.0)/(2.0.*sqrt(-gtt[iy,iz])))             
        end
    end
    return nrg
end

function energyNLC(Q,Qdx,Qdy,Qdz,Qdot)
    nrg=0.0;
    @inbounds for iy=2:yvaltest
        @inbounds for iz=zvaltest1:zvaltest2
                fact=dy*dz*1.0/sqrt(1.0 - (tan((pi*ys[iy])/2.0)^2.0 + tan((pi*zs[iz])/2.0)^2.0)/(a + b*(tan((pi*ys[iy])/2.0)^2.0 + tan((pi*zs[iz])/2.0)^2.0)^2.0))*(2.0*pi*tan(pi/2.0*ys[iy]))*(pi^2.0/4.0*sec(pi*ys[iy]/2.0)^2.0*sec(pi*zs[iz]/2.0)^2.0)
                
                nrg+=fact*((2.0*Q[iy,iz]^(1.0 + 3.0) + (1.0 + 3.0)*(gzz[iy,iz]*Qdz[iy,iz]^2.0 + Qdy[iy,iz]*(2.0*gyz[iy,iz]*Qdz[iy,iz] + gyy[iy,iz]*Qdy[iy,iz]) + 2.0*(gxz[iy,iz]*Qdz[iy,iz] + gxy[iy,iz]*Qdy[iy,iz])*Qdx[iy,iz] + gxx[iy,iz]*Qdx[iy,iz]^2.0 - gtt[iy,iz]*Qdot[iy,iz]^2.0))/(2.0.*(1.0 + 3.0)*sqrt(-gtt[iy,iz])))          
        end
    end
    return nrg
end
