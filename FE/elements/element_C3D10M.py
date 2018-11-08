from elements.element_base import Element_Base
import numpy as np


class Element_C3D10M(Element_Base):
    
    def __init__(self):
        # coordinates - can be reference or current
        # 11th coord is the auxialliary 
        Element_Base.__init__(self,3,10,3)
        
        self.nsdm   = self.ndim*self.ndof             # stress dimension = 3*3 = 9
        self.nvoigt = int(self.ndim*(self.ndim+1)/2)  # number of stress terms in voigt notation = 6 
        self.nquad  = 5
        self.nsub = [
                [0,4,6,7],
                [1,5,4,8],
                [2,6,5,9],
                [3,8,7,9],
                [4,8,5,10],
                [5,8,9,10],
                [9,8,7,10],
                [7,8,4,10],
                [4,5,6,10],
                [5,9,6,10],
                [9,7,6,10],
                [7,4,6,10]
            ]
        
        # ============================
        # coords of integration points
        # sg is 4xnquad
        # see Table II in P. THOUTIREDDY ET AL
        # ============================
        self.sg = [
                [0.25,0.25,0.25,0.25],
                [0.5,0.1666666666666667,0.1666666666666667,0.1666666666666667],
                [0.1666666666666667,0.5,0.1666666666666667,0.1666666666666667],
                [0.1666666666666667,0.1666666666666667,0.5,0.1666666666666667],
                [0.1666666666666667,0.1666666666666667,0.1666666666666667,0.5]
            ]
        
        # convert list to numpy array
        self.sg = np.asarray(self.sg) # 5 x 4
        self.sg = self.sg.transpose() # 4 x 5
        # sg(1,1) = 0.25
        # sg(2,1) = 0.25
        # sg(3,1) = 0.25
        # sg(4,1) = 0.25 = 1 - s1 - s2 - s3
        #
        # sg(1,2) = 0.5
        # sg(2,2) = 1./6.
        # sg(3,2) = 1./6.
        # sg(4,2) = 1./6. = 1 - s1 - s2 - s3
        #
        # sg(1,3) = 1./6.
        # sg(2,3) = 0.5
        # sg(3,3) = 1./6.
        # sg(4,3) = 1./6. = 1 - s1 - s2 - s3
        #
        # sg(1,4) = 1./6.
        # sg(2,4) = 1./6.
        # sg(3,4) = 0.5
        # sg(4,4) = 1./6. = 1 - s1 - s2 - s3
        #
        # sg(1,5) = 1./6.
        # sg(2,5) = 1./6.
        # sg(3,5) = 1./6.
        # sg(4,5) = 0.5   = 1 - s1 - s2 - s3
        
        
        # ============================
        # weights
        # weights is 12xnquad
        # below are the actually the multipliers of the subtet volume
        # see Table II in P. THOUTIREDDY ET AL        
        # ============================
        self.weights = [
                # V1 V2 V3 V4  V5 V6  V7  V8  V9 V10 V11 V12
                [ 0, 0, 0, 0, 0.5, 1, 0.5, 1, 1, 0.5, 1, 0.5], # quad1
                [ 1, 0, 0, 0,   0, 0,   0, 0, 0,   0, 0, 0.5], # quad2
                [ 0, 1, 0, 0, 0.5, 0,   0, 0, 0,   0, 0,   0], # quad3
                [ 0, 0, 1, 0,   0, 0,   0, 0, 0, 0.5, 0,   0], # quad4
                [ 0, 0, 0, 1,   0, 0, 0.5, 0, 0,   0, 0,   0]  # quad5
            ]
        # convert list to numpy array
        self.weights = np.asarray(self.weights)  # 5 x 12
        self.weights = self.weights.transpose()  # 12 x 5
        
        self.xcurrent = np.zeros((3,11))
        self.xref     = np.zeros((3,11))
        
        
        # ------------------------------------------------------------------
        # gradient matrices, tangent moduli and stress vector (voigt-style)
	    # ------------------------------------------------------------------
        self.Bmatl = np.zeros((self.nquad,self.nvoigt,self.nnodes*self.ndof)) # nvoigt = 3 for 2d (plane-strain/stress), 6 for 3d
        self.Gmatl = np.zeros((self.nquad,self.nsdm,self.nnodes*self.ndof))
        self.Amatl = np.zeros((self.nquad,self.nsdm,self.nsdm))               # note the size of Dmat is for general unsymmetric moduli, based on D.22 pg 763 in Peric complas      
        self.Svecl = np.zeros((self.nquad,self.nvoigt))                       # note: symmetric cauchy stress
	



        
    def compute_shape(self):
        
        
        # compute the derivatives
        for lquad in range(0,self.nquad): #{
            
            a = np.zeros((4,4)) 
        
            # shape functions
            # shp(0,:) = derivs
            # shp(1,:) = derivs
            # shp(2,:) = derivs
            # shp(3,:) = actual shape functions
            shp = np.zeros((4,10))
            
            bx  = np.zeros((4,10))
            by  = np.zeros((4,10))
            bz  = np.zeros((4,10))
            xsj = np.zeros(self.nquad)
            
        
            shp[3,0] = self.sg[0,lquad]*(2.*self.sg[0,lquad]-1.)
            shp[3,1] = self.sg[1,lquad]*(2.*self.sg[1,lquad]-1.)
            shp[3,2] = self.sg[2,lquad]*(2.*self.sg[2,lquad]-1.)
            shp[3,3] = self.sg[3,lquad]*(2.*self.sg[3,lquad]-1.)
            shp[3,4] = 4.*self.sg[0,lquad]*self.sg[1,lquad]
            shp[3,5] = 4.*self.sg[1,lquad]*self.sg[2,lquad]
            shp[3,6] = 4.*self.sg[2,lquad]*self.sg[0,lquad]
            shp[3,7] = 4.*self.sg[0,lquad]*self.sg[3,lquad]
            shp[3,8] = 4.*self.sg[1,lquad]*self.sg[3,lquad]
            shp[3,9] = 4.*self.sg[2,lquad]*self.sg[3,lquad]
        
            
            # compute the auxilliary node coordinate
            self.xref[:,10]=(self.xref[:,4]+self.xref[:,5]+self.xref[:,6]+\
                            self.xref[:,7]+self.xref[:,8]+self.xref[:,9])/6.
        
            # ===========================================
            # compute derivatives of shape functions
            # ie we want to compute the term \bar{L}_{aJ}
            # see eq (18) in P. THOUTIREDDY ET AL
            # ===========================================
            # -------------------------------------------------------------------
            # loop through the 12 sub-tetras
            xj=np.zeros(12) # array of 6V's of sub-tetras
            for i in range(0,12): #{
                # 
                # compute the derivatives of shape functions of the individual sub-tetras
                #
                # notes: 
                #  1) We use the linear shape functions:
                #      N1 = s1
                #      N2 = s2
                #      N3 = s3
                #      N4 = 1-s1-s2-s3
                #    Therefore,
                #
                #     dN/ds_1 = [1,0,0,-1]
                #     dN/ds_2 = [0,1,0,-1]
                #     dN/ds_3 = [0,0,1,-1]
                #
                #     ds(3,4) == dN/ds
                
                #
                #  2) x = N_1*x_1 + N_2*x_2 + N_3*x_3 + N_4*x_4
                #     dX/ds_l = dN_1/ds_l*x_1 + dN_2/ds_l*x_2 + dN_3/ds_l*x_3 + dN_4/ds_l*x_4
                #             = [dN_1/ds_l + dN_2/ds_l + dN_3/ds_l  dN_4/ds_l] * [x1 x2 x3 x4]^T
                #             
                #     Therefore, a(3,3) == dX/ds where:
                #
                #     dX/ds_1 = x1-x4
                #     dX/ds_2 = x2-x4
                #     dX/ds_3 = x3-x4
                #
                #     dY/ds_1 = x1-x4
                #     dY/ds_2 = x2-x4
                #     dY/ds_3 = x3-x4
                #
                #     dZ/ds_1 = x1-x4
                #     dZ/ds_2 = x2-x4
                #     dZ/ds_3 = x3-x4
                #
                #
                #  3) dN/dX = (ds/dX)*(dN/ds)
                #           = (dX/ds)^{-1}*(dN/ds)
                #           = (3x3)*(3x4)
                #           = (3x4)
                #     or:
                #      [dx]   [dN_1/dX dN_2/dX dN_3/dX dN_4/dX]
                #      [dy] = [dN_1/dY dN_2/dY dN_3/dY dN_4/dY]
                #      [dz]   [dN_1/dZ dN_2/dZ dN_3/dZ dN_4/dZ]
                #
                
                dx = np.zeros((11,12))
                dy = np.zeros((11,12))
                dz = np.zeros((11,12))
                
                
                # nodes of subtet
                n1 = self.nsub[i][0]
                n2 = self.nsub[i][1]
                n3 = self.nsub[i][2]
                n4 = self.nsub[i][3]
                
                # coords of subtet
                x1 = self.xref[0,n1]
                y1 = self.xref[1,n1]
                z1 = self.xref[2,n1]
                
                x2 = self.xref[0,n2]
                y2 = self.xref[1,n2]
                z2 = self.xref[2,n2]
                
                x3 = self.xref[0,n3]
                y3 = self.xref[1,n3]
                z3 = self.xref[2,n3]
                
                x4 = self.xref[0,n4]
                y4 = self.xref[1,n4]
                z4 = self.xref[2,n4]
                
                dNds = np.array([
                                 [  1.0,  0.0,  0.0,  -1.0], 
                                 [  0.0,  1.0,  0.0,  -1.0], 
                                 [  0.0,  0.0,  1.0,  -1.0]
                                ]) # 3x4
                
                temp = np.array([
                                 [1.0,1.0,1.0,1.0], 
                                 [ x1, x2, x3, x4],
                                 [ y1, y2, y3, y4],
                                 [ z1, z2, z3, z4]
                                ])
                
                xj[i] = np.linalg.det(temp) # this is 6*V_subtet
                
                dxds = np.array([
                     [x1-x4,  x2-x4,  x3-x4], 
                     [y1-y4,  y2-y4,  y3-y4], 
                     [z1-z4,  z2-z4,  z3-z4]
                     ])
    
            
                dsdx = np.linalg.inv(dxds)
                
                dNdx   = np.matmul(dsdx,dNds) # 3x4
                
                dx[n1,i] = dNdx[0,0]
                dy[n1,i] = dNdx[1,0]
                dz[n1,i] = dNdx[2,0]
                
                dx[n2,i] = dNdx[0,1]
                dy[n2,i] = dNdx[1,1]
                dz[n2,i] = dNdx[2,1]
                
                dx[n3,i] = dNdx[0,2]
                dy[n3,i] = dNdx[1,2]
                dz[n3,i] = dNdx[2,2]
                
                dx[n4,i] = dNdx[0,3]
                dy[n4,i] = dNdx[1,3]
                dz[n4,i] = dNdx[2,3]
                
            #} end loop of sub-tetras
            # -------------------------------------------------------------------
            
            
            # -------------------------------------------------------------------
            # loop over quadrature points of the composite tetra
            xsj = np.zeros(self.nquad)
            a   = np.zeros((4,4))
            for l in range(0,self.nquad): #{
                
                
                lam0 = 1.0
                lam1 = self.sg[0,l]-0.25
                lam2 = self.sg[1,l]-0.25
                lam3 = self.sg[2,l]-0.25
                
                
                for j in range(0,12):
                    xsj[l] +=  np.absolute( xj[j] )*self.weights[j,l]/6.
                
                
                # compute the term given by eq (19)
                for j in range(0,12):
                    
                    
                    a[0,0] += xj[j]*self.weights[j,l]*lam0
                    a[0,1] += xj[j]*self.weights[j,l]*lam1
                    a[0,2] += xj[j]*self.weights[j,l]*lam2
                    a[0,3] += xj[j]*self.weights[j,l]*lam3  
                    
                    a[1,1] += xj[j]*self.weights[j,l]*lam1*lam1
                    a[1,2] += xj[j]*self.weights[j,l]*lam2*lam1
                    a[1,3] += xj[j]*self.weights[j,l]*lam3*lam1  
                    
                    a[2,2] += xj[j]*self.weights[j,l]*lam2*lam2
                    a[2,3] += xj[j]*self.weights[j,l]*lam3*lam2
                    
                    a[3,3] += xj[j]*self.weights[j,l]*lam3*lam3
                    
                for j in range(0,12):
                    
                    for i in range(0,10):
                        bx[0,i] +=  xj[j]*self.weights[j,l]*dx[i,j]*lam0
                        bx[1,i] +=  xj[j]*self.weights[j,l]*dx[i,j]*lam1
                        bx[2,i] +=  xj[j]*self.weights[j,l]*dx[i,j]*lam2
                        bx[3,i] +=  xj[j]*self.weights[j,l]*dx[i,j]*lam3
                        
                        by[0,i] +=  xj[j]*self.weights[j,l]*dx[i,j]*lam0
                        by[1,i] +=  xj[j]*self.weights[j,l]*dx[i,j]*lam1
                        by[2,i] +=  xj[j]*self.weights[j,l]*dx[i,j]*lam2
                        by[3,i] +=  xj[j]*self.weights[j,l]*dx[i,j]*lam3
                        
                        bz[0,i] +=  xj[j]*self.weights[j,l]*dx[i,j]*lam0
                        bz[1,i] +=  xj[j]*self.weights[j,l]*dx[i,j]*lam1
                        bz[2,i] +=  xj[j]*self.weights[j,l]*dx[i,j]*lam2
                        bz[3,i] +=  xj[j]*self.weights[j,l]*dx[i,j]*lam3
                        
                    for i in range(4,10):
                        
                        bx[0,i] +=  xj[j]*self.weights[j,l]*dx[10,j]*lam0/6.0
                        bx[1,i] +=  xj[j]*self.weights[j,l]*dx[10,j]*lam1/6.0
                        bx[2,i] +=  xj[j]*self.weights[j,l]*dx[10,j]*lam2/6.0
                        bx[3,i] +=  xj[j]*self.weights[j,l]*dx[10,j]*lam3/6.0
                        
                        by[0,i] +=  xj[j]*self.weights[j,l]*dx[10,j]*lam0/6.0
                        by[1,i] +=  xj[j]*self.weights[j,l]*dx[10,j]*lam1/6.0
                        by[2,i] +=  xj[j]*self.weights[j,l]*dx[10,j]*lam2/6.0
                        by[3,i] +=  xj[j]*self.weights[j,l]*dx[10,j]*lam3/6.0
                        
                        bz[0,i] +=  xj[j]*self.weights[j,l]*dx[10,j]*lam0/6.0
                        bz[1,i] +=  xj[j]*self.weights[j,l]*dx[10,j]*lam1/6.0
                        bz[2,i] +=  xj[j]*self.weights[j,l]*dx[10,j]*lam2/6.0
                        bz[3,i] +=  xj[j]*self.weights[j,l]*dx[10,j]*lam3/6.0
                        
            #} end loop of quadrature points of the composite tetra
            # -------------------------------------------------------------------
            a[1,0]=a[0,1]
            a[2,0]=a[0,2]
            a[2,1]=a[1,2]
            
            a[3,0]=a[0,3]
            a[3,1]=a[1,3]
            a[3,2]=a[2,3]
    
            ai = np.linalg.inv(a)
            
            a0 = np.zeros((3,10))
            a1 = np.zeros((3,10))
            a2 = np.zeros((3,10))
            a3 = np.zeros((3,10))
            
            
            for i in range(0,10):
                for j in range(0,4):
                    
                    a0[0,i] += ai[0,j]*bx[j,i]
                    a0[1,i] += ai[0,j]*by[j,i]
                    a0[2,i] += ai[0,j]*bz[j,i]
                    
                    a1[0,i] += ai[1,j]*bx[j,i]
                    a1[1,i] += ai[1,j]*by[j,i]
                    a1[2,i] += ai[1,j]*bz[j,i]
                    
                    a2[0,i] += ai[2,j]*bx[j,i]
                    a2[1,i] += ai[2,j]*by[j,i]
                    a2[2,i] += ai[2,j]*bz[j,i]
                    
                    a3[0,i] += ai[3,j]*bx[j,i]
                    a3[1,i] += ai[3,j]*by[j,i]
                    a3[2,i] += ai[3,j]*bz[j,i]
                    
            
            lam0 = 1.0
            lam1 = self.sg[0,lquad]-0.25
            lam2 = self.sg[1,lquad]-0.25
            lam3 = self.sg[2,lquad]-0.25
                
                
            for i in range(0,10):
                for j in range(0,3):
                    shp[j,i]=a0[j,i]*lam0+a1[j,i]*lam1+a2[j,i]*lam2+a3[j,i]*lam3
    
            
            # assign to matrices
            for a in range(0,self.nnodes):
                
                i0 = a*self.ndof
                i1 = i0+1
                i2 = i0+2
                # -----------------------------
                # Displacement-Strain operator
                # -----------------------------
                self.Bmatl[lquad,0,i0] = shp[0,a] # d/dx -> 11
                self.Bmatl[lquad,1,i1] = shp[1,a] # d/dy -> 22
                self.Bmatl[lquad,2,i2] = shp[2,a] # d/dz -> 33
                
                self.Bmatl[lquad,3,i0] = shp[1,a] # d/dy -> 12
                self.Bmatl[lquad,3,i1] = shp[0,a] # d/dx -> 21
                
                self.Bmatl[lquad,4,i0] = shp[2,a] # d/dz -> 13
                self.Bmatl[lquad,4,i2] = shp[0,a] # d/dx -> 31
                
                self.Bmatl[lquad,3,i1] = shp[2,a] # d/dz -> 23
                self.Bmatl[lquad,3,i2] = shp[1,a] # d/dy -> 32
                
                # -----------------------------
                # Gradient operator
                # -----------------------------
                self.Gmatl[lquad,0,i0] = shp[0,a] # d/dx -> 11
                self.Gmatl[lquad,1,i1] = shp[0,a] # d/dx -> 11
                self.Gmatl[lquad,2,i2] = shp[0,a] # d/dx -> 11
                
                self.Gmatl[lquad,3,i0] = shp[1,a] # d/dy -> 22
                self.Gmatl[lquad,4,i1] = shp[1,a] # d/dy -> 22
                self.Gmatl[lquad,5,i2] = shp[1,a] # d/dy -> 22
                
                self.Gmatl[lquad,6,i0] = shp[2,a] # d/dz -> 33
                self.Gmatl[lquad,7,i1] = shp[2,a] # d/dz -> 33
                self.Gmatl[lquad,8,i2] = shp[2,a] # d/dz -> 33
                
                

            
        #} end loop of gauss points
        
        