#RK steps, dissipations and data storage
#Author: Alejandro Cardenas-Avendano

function simulation!(steps,dtR,phi_Mn1,pi_Mn1,phi_Mn2,pi_Mn2,phidx_Mn,pidx_Mn,phidy_Mn,pidy_Mn,phidz_Mn,pidz_Mn,phidot_Mn,pidot_Mn,phidxdx_Mn,phidxdy_Mn,phidxdz_Mn,phidydy_Mn,phidydz_Mn,phidzdz_Mn)
    aux=1
    println("Number of steps:",steps);

    for t in 1:steps

        #Compute the RHS
        RHS_Eq!(phi_Mn1,pi_Mn1,phidx_Mn,pidx_Mn,phidy_Mn,pidy_Mn,phidz_Mn,pidz_Mn,phidot_Mn,pidot_Mn,phidxdx_Mn,phidxdy_Mn,phidxdz_Mn,phidydy_Mn,phidydz_Mn,phidzdz_Mn)

        #Save snapshot
        if (ts_aux[t]%dtR)==0

            if (ts_aux[t]%dtRFS)==0

                println("Saving Snapshot at ",ts[t],"/",Tf);

                @inbounds for iz=2:(Nz0-1)
                    phidot_M[1,iz]=-pi_Mn1[1,iz]/gtt[1,iz]
                end
    
                spatial!(phidot_M,phidotdy_M,phidotdz_M)

                energy=energyNLC(phi_Mn1,phidx_M,phidy_M,phidz_M,phidot_M)

                fileData=location*"Wave_$(filename)_$(Am)_$(ts[t])_$(dy)_$(lambda)_$(epsKO).h5"
                # #=
                if isfile(fileData)
                    rm(fileData)
                    println("File Data Overwritten")
                end
                
                #These two are mandatory to store! Used to continue an evolution
                h5write(fileData, "phi",phi_Mn1[:,:])
                h5write(fileData, "pi",pi_Mn1[:,:])
                
                h5write(fileData, "energy",energy)
               
                println("");
                println("Energy ",energy);
                println("");
            end

            if isfile(fileMax)
                rm(fileMax)
                println("File Maxs Overwritten")
            end
            
            h5write(fileMax, "ts",ts[:])
            h5write(fileMax, "pidotmaxs",pidotmaxval[:])
            h5write(fileMax, "phidthdthmaxs",phidthdthmaxval[:])

            aux+=1;
            
        end

        pidotmaxval[t]=maximum(abs.(pidot_Mn[1:yvaltest,zvaltest1:zvaltest2]))


        @inbounds Threads.@threads for iy=2:(Ny0-1)
            @inbounds for iz=2:(Nz0-1)
                    phik1_M[iy,iz] = (((-π*(1.0+cos(π*ys[iy]))*sqrt(tan(π*ys[iy]/2.0)^2.0+tan(π*zs[iz]/2.0)^2.0))*phidy_Mn[iy,iz]/sqrt(1.0+cot(π*ys[iy]/2.0)^2.0*tan(π*zs[iz]/2.0)^2.0))+4.0*tan(π*ys[iy]/2.0)^2.0*phidxdx_Mn[iy,iz])/π^2.0;
                    pik1_M[iy,iz] = ((pi*(-2.0 + (-1.0 + cos(pi*ys[iy]))*cos(pi*zs[iz]))*sin(pi*zs[iz])*phidz_Mn[iy,iz])/(1.0 + cos(pi*ys[iy])) + (2.0*(4.0*cos((pi*zs[iz])/2.0)^6.0*tan((pi*ys[iy])/2.0)^2.0*phidzdz_Mn[iy,iz] + (sin((pi*ys[iy])/2.0)^2.0*sqrt(1.0 + cot((pi*ys[iy])/2.0)^2.0*tan((pi*zs[iz])/2.0)^2.0)*(pi*(-2.0 + cos(pi*ys[iy])*(-1.0 + cos(pi*zs[iz])))*phidy_Mn[iy,iz] - (2.0*sin(pi*zs[iz]) + sin(2.0*pi*zs[iz]))*phidydz_Mn[iy,iz]))/sqrt(tan((pi*ys[iy])/2.0)^2.0 + tan((pi*zs[iz])/2.0)^2.0)+ 4.0*cos((pi*ys[iy])/2.0)^4.0*sin((pi*zs[iz])/2.0)^2.0*phidydy_Mn[iy,iz]))/(1.0 + cos(pi*zs[iz])))/pi^2.0;
            end
        end

        phidthdthmaxval[t]=maximum(abs.(pik1_M[1:yvaltest,zvaltest1:zvaltest2]))

        #Fourth order Runge-Kutta for the time integration

        #Calculate k1 term in RK4 formula
        @inbounds Threads.@threads for iy=2:(Ny0-1)
            @inbounds for iz=2:(Nz0-1)
                    phik1_M[iy,iz] = phidot_Mn[iy,iz];
                    pik1_M[iy,iz]  = pidot_Mn[iy,iz];
            end
        end

        #Calculate k2 term in RK4 formula

        @inbounds Threads.@threads for iy=2:(Ny0-1)
            @inbounds for iz=2:(Nz0-1)
                    phi_Mn2[iy,iz] = phi_Mn1[iy,iz] + 0.5 * dt * phik1_M[iy,iz];
                    pi_Mn2[iy,iz]  = pi_Mn1[iy,iz]  + 0.5 * dt * pik1_M[iy,iz];
            end
        end

        
        @inbounds Threads.@threads for iz=2:(Nz0-1)
            phi_Mn2[1,iz] = fw_diffC4(phi_Mn2,1,iz,fw1yC4);
            #pi_Mn2[1,iz]  = fw_diffC4(pi_Mn2,1,iz,fw1yC4);
        end

        RHS_Eq!(phi_Mn2,pi_Mn2,phidx_Mn,pidx_Mn,phidy_Mn,pidy_Mn,phidz_Mn,pidz_Mn,phidot_Mn,pidot_Mn,phidxdx_Mn,phidxdy_Mn,phidxdz_Mn,phidydy_Mn,phidydz_Mn,phidzdz_Mn)

        @inbounds Threads.@threads for iy=2:(Ny0-1)
            @inbounds for iz=2:(Nz0-1)
                    phik2_M[iy,iz] =  phidot_Mn[iy,iz];
                    pik2_M[iy,iz]  =  pidot_Mn[iy,iz];
            end
        end

        #Calculate k3 term in RK4 formula

        @inbounds Threads.@threads for iy=2:(Ny0-1)
            @inbounds for iz=2:(Nz0-1)
                    phi_Mn2[iy,iz] = phi_Mn1[iy,iz] + 0.5 * dt * phik2_M[iy,iz];
                    pi_Mn2[iy,iz]  = pi_Mn1[iy,iz]  + 0.5 * dt * pik2_M[iy,iz];
            end
        end

        @inbounds Threads.@threads for iz=2:(Nz0-1)
            phi_Mn2[1,iz] = fw_diffC4(phi_Mn2,1,iz,fw1yC4);
            #pi_Mn2[1,iz]  = fw_diffC4(pi_Mn2,1,iz,fw1yC4);
        end

        RHS_Eq!(phi_Mn2,pi_Mn2,phidx_Mn,pidx_Mn,phidy_Mn,pidy_Mn,phidz_Mn,pidz_Mn,phidot_Mn,pidot_Mn,phidxdx_Mn,phidxdy_Mn,phidxdz_Mn,phidydy_Mn,phidydz_Mn,phidzdz_Mn)
        
        @inbounds Threads.@threads for iy=2:(Ny0-1)
            @inbounds for iz=2:(Nz0-1)
                    phik3_M[iy,iz] = phidot_Mn[iy,iz];
                    pik3_M[iy,iz]  = pidot_Mn[iy,iz];
            end
        end

        #Calculate k4 term in RK4 formula

        @inbounds Threads.@threads for iy=2:(Ny0-1)
            @inbounds for iz=2:(Nz0-1)
                    phi_Mn2[iy,iz] = phi_Mn1[iy,iz] + dt * phik3_M[iy,iz];
                    pi_Mn2[iy,iz]  = pi_Mn1[iy,iz]  + dt * pik3_M[iy,iz];
            end
        end

        @inbounds Threads.@threads for iz=2:(Nz0-1)
            phi_Mn2[1,iz] = fw_diffC4(phi_Mn2,1,iz,fw1yC4);
        end

        RHS_Eq!(phi_Mn2,pi_Mn2,phidx_Mn,pidx_Mn,phidy_Mn,pidy_Mn,phidz_Mn,pidz_Mn,phidot_Mn,pidot_Mn,phidxdx_Mn,phidxdy_Mn,phidxdz_Mn,phidydy_Mn,phidydz_Mn,phidzdz_Mn)

        @inbounds Threads.@threads for iy=2:(Ny0-1)
            @inbounds for iz=2:(Nz0-1)
                    phik4_M[iy,iz] = phidot_Mn[iy,iz];
                    pik4_M[iy,iz]  = pidot_Mn[iy,iz];
            end
        end

        #Time evolution RK4 formula

        @inbounds Threads.@threads for iy=2:(Ny0-1)
            @inbounds for iz=2:(Nz0-1)
                    phi_Mn1[iy,iz] = phi_Mn1[iy,iz] + dt * (phik1_M[iy,iz] + 2.0 * phik2_M[iy,iz] + 2.0 * phik3_M[iy,iz] + phik4_M[iy,iz]) / 6.0;
                    pi_Mn1[iy,iz]  = pi_Mn1[iy,iz]  + dt * (pik1_M[iy,iz]  + 2.0 * pik2_M[iy,iz]  + 2.0 * pik3_M[iy,iz]  + pik4_M[iy,iz])  / 6.0;
            end
        end

        @inbounds Threads.@threads for iz=2:(Nz0-1)
            phi_Mn1[1,iz] = fw_diffC4(phi_Mn1,1,iz,fw1yC4);
            pi_Mn1[1,iz]  = fw_diffC4(pi_Mn1,1,iz,fw1yC4);
            phidot_M[1,iz]=-pi_Mn1[1,iz]/gtt[1,iz]
        end

        #########################
        #Kreiss-Oliger dissipation 
        #########################

        #Zero phi_Mn2 and  pi_Mn2 as auxiliar variables to compute dissipation
        #This may not be required
        
        # #=
        @inbounds Threads.@threads for iy=1:(Ny0)
            @inbounds for iz=1:(Nz0)   
                phi_Mn2[iy,iz]=0.0
                pi_Mn2[iy,iz] =0.0 
            end
        end
        # =#
        
        #Center of the grid
        # #=
        @inbounds Threads.@threads for iy=4:(Ny0-3)
            @inbounds for iz=4:(Nz0-3)  
                phi_Mn2[iy,iz]=cd_diss6(phi_Mn1,iy,iz,cd6y)+cd_diss6(phi_Mn1,iy,iz,cd6z)
                pi_Mn2[iy,iz]=cd_diss6(pi_Mn1,iy,iz,cd6y)+cd_diss6(pi_Mn1,iy,iz,cd6z)
            end
        end
        # =#

        # #=
        @inbounds Threads.@threads for iz=4:(Nz0-3)

            phi_Mn2[3,iz]=cd_diss6(phi_Mn1,3,iz,cd6yL1)+cd_diss6(phi_Mn1,3,iz,cd6z)
            pi_Mn2[3,iz]=cd_diss6(pi_Mn1,3,iz,cd6yL1)+cd_diss6(pi_Mn1,3,iz,cd6z)

            phi_Mn2[2,iz]=cd_diss6(phi_Mn1,2,iz,cd6yL2)+cd_diss6(phi_Mn1,2,iz,cd6z)
            pi_Mn2[2,iz]=cd_diss6(pi_Mn1,2,iz,cd6yL2)+cd_diss6(pi_Mn1,2,iz,cd6z)

            phi_Mn2[1,iz]=cd_diss6(phi_Mn1,1,iz,cd6yL3)+cd_diss6(phi_Mn1,1,iz,cd6z)
            pi_Mn2[1,iz]=cd_diss6(pi_Mn1,1,iz,cd6yL3)+cd_diss6(pi_Mn1,1,iz,cd6z)

            phi_Mn2[Ny0-1,iz]=cd_diss6(phi_Mn1,Ny0-1,iz,cd6yR2)+cd_diss6(phi_Mn1,Ny0-1,iz,cd6z)
            pi_Mn2[Ny0-1,iz]=cd_diss6(pi_Mn1,Ny0-1,iz,cd6yR2)+cd_diss6(pi_Mn1,Ny0-1,iz,cd6z)

            phi_Mn2[Ny0-2,iz]=cd_diss6(phi_Mn1,Ny0-2,iz,cd6yR1)+cd_diss6(phi_Mn1,Ny0-2,iz,cd6z)
            pi_Mn2[Ny0-2,iz]=cd_diss6(pi_Mn1,Ny0-2,iz,cd6yR1)+cd_diss6(pi_Mn1,Ny0-2,iz,cd6z)

        end
        # =#

        phi_Mn2[1,2]=cd_diss6(phi_Mn1,1,2,cd6yL3)+cd_diss6(phi_Mn1,1,2,cd6zL2)
        pi_Mn2[1,2]=cd_diss6(pi_Mn1,1,2,cd6yL3)+cd_diss6(pi_Mn1,1,2,cd6zL2)

        phi_Mn2[2,3]=cd_diss6(phi_Mn1,2,3,cd6yL2)+cd_diss6(phi_Mn1,2,3,cd6zL1)
        pi_Mn2[2,3]=cd_diss6(pi_Mn1,2,3,cd6yL2)+cd_diss6(pi_Mn1,2,3,cd6zL1)

        phi_Mn2[1,Nz0-1]=cd_diss6(phi_Mn1,1,Nz0-1,cd6yL3)+cd_diss6(phi_Mn1,1,Nz0-1,cd6zR2)
        pi_Mn2[1,Nz0-1]=cd_diss6(pi_Mn1,1,Nz0-1,cd6yL3)+cd_diss6(pi_Mn1,1,Nz0-1,cd6zR2)

        phi_Mn2[2,Nz0-2]=cd_diss6(phi_Mn1,2,Nz0-2,cd6yL2)+cd_diss6(phi_Mn1,2,Nz0-2,cd6zR1)
        pi_Mn2[2,Nz0-2]=cd_diss6(pi_Mn1,2,Nz0-2,cd6yL2)+cd_diss6(pi_Mn1,2,Nz0-2,cd6zR1)
        
        @inbounds Threads.@threads for iy=4:(Ny0-4)

            phi_Mn2[iy,2]=cd_diss6(phi_Mn1,iy,2,cd6y)+cd_diss6(phi_Mn1,iy,2,cd6zL2)
            pi_Mn2[iy,2]=cd_diss6(pi_Mn1,iy,2,cd6y)+cd_diss6(pi_Mn1,iy,2,cd6zL2)

            phi_Mn2[iy,3]=cd_diss6(phi_Mn1,iy,3,cd6y)+cd_diss6(phi_Mn1,iy,3,cd6zL1)
            pi_Mn2[iy,3]=cd_diss6(pi_Mn1,iy,3,cd6y)+cd_diss6(pi_Mn1,iy,3,cd6zL1)

            phi_Mn2[iy,Nz0-1]=cd_diss6(phi_Mn1,iy,Nz0-1,cd6y)+cd_diss6(phi_Mn1,iy,Nz0-1,cd6zR2)
            pi_Mn2[iy,Nz0-1]=cd_diss6(pi_Mn1,iy,Nz0-1,cd6y)+cd_diss6(pi_Mn1,iy,Nz0-1,cd6zR2)

            phi_Mn2[iy,Nz0-2]=cd_diss6(phi_Mn1,iy,Nz0-2,cd6y)+cd_diss6(phi_Mn1,iy,Nz0-2,cd6zR1)
            pi_Mn2[iy,Nz0-2]=cd_diss6(pi_Mn1,iy,Nz0-2,cd6y)+cd_diss6(pi_Mn1,iy,Nz0-2,cd6zR1)

        end

        #Application of the filter
        @inbounds Threads.@threads for iy=2:Ny0-1
            @inbounds for iz=2:Nz0-1
                phi_Mn1[iy,iz]+=phi_Mn2[iy,iz]
                pi_Mn1[iy,iz] +=pi_Mn2[iy,iz]   
            end
        end

        if any(x->x>1e4, phi_Mn1)
            println("No bueno! ")
            break;
        end

    end
end
