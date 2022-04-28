#include "mex.h" 
#include <stdio.h>
#include <math.h>
#include <string.h>

#define TWOPI   6.28318530717959
#define GAMMA   TWOPI
#define DEBUG

// CTR:
#define DEBUG_printf( ... ) \
if (debugflag) mexPrintf(__VA_ARGS__); \
else (void)0
// (The (void)0 gives an error if DEBUG_printf is missing a terminating ;.
static bool debugflag = false;
// End CTR.


/* Multiply 3x3 matrix by 3x1 vector. */
void multmatvec(double *mat, double *vec, double *matvec)
{
    *matvec++ = mat[0]*vec[0] + mat[3]*vec[1] + mat[6]*vec[2];
    *matvec++ = mat[1]*vec[0] + mat[4]*vec[1] + mat[7]*vec[2];
    *matvec++ = mat[2]*vec[0] + mat[5]*vec[1] + mat[8]*vec[2];
}


/* Add two 3x1 Vectors */
void addvecs(double *vec1, double *vec2, double *vecsum)
{
    *vecsum++ = *vec1++ + *vec2++;
    *vecsum++ = *vec1++ + *vec2++;
    *vecsum++ = *vec1++ + *vec2++;
}


/* ======== Adjoint of a 3x3 matrix ========= */
void adjmat(double *mat, double *adj)
{
    *adj++ = (mat[4]*mat[8]-mat[7]*mat[5]);
    *adj++ =-(mat[1]*mat[8]-mat[7]*mat[2]);
    *adj++ = (mat[1]*mat[5]-mat[4]*mat[2]);
    *adj++ =-(mat[3]*mat[8]-mat[6]*mat[5]);
    *adj++ = (mat[0]*mat[8]-mat[6]*mat[2]);
    *adj++ =-(mat[0]*mat[5]-mat[3]*mat[2]);
    *adj++ = (mat[3]*mat[7]-mat[6]*mat[4]);
    *adj++ =-(mat[0]*mat[7]-mat[6]*mat[1]);
    *adj++ = (mat[0]*mat[4]-mat[3]*mat[1]);
}


/* ====== Set a 3x3 matrix to all zeros	======= */
void zeromat(double *mat)
{
    *mat++=0;
    *mat++=0;
    *mat++=0;
    *mat++=0;
    *mat++=0;
    *mat++=0;
    *mat++=0;
    *mat++=0;
    *mat++=0;
}


/* ======== Return 3x3 Identity Matrix  ========= */
void eyemat(double *mat)
{
    zeromat(mat);
    mat[0]=1;
    mat[4]=1;
    mat[8]=1;
}


/* ======== Determinant of a 3x3 matrix ======== */
double detmat(double *mat)
{
    double det;
    
    det = mat[0]*mat[4]*mat[8];
    det+= mat[3]*mat[7]*mat[2];
    det+= mat[6]*mat[1]*mat[5];
    det-= mat[0]*mat[7]*mat[5];
    det-= mat[3]*mat[1]*mat[8];
    det-= mat[6]*mat[4]*mat[2];
    
    return det;
}


/* ======== multiply a matrix by a scalar ========= */
void scalemat(double *mat, double scalar)
{
    *mat++ *= scalar;
    *mat++ *= scalar;
    *mat++ *= scalar;
    *mat++ *= scalar;
    *mat++ *= scalar;
    *mat++ *= scalar;
    *mat++ *= scalar;
    *mat++ *= scalar;
    *mat++ *= scalar;
}


/* ======== Inverse of a 3x3 matrix ========= */
/*	DO NOT MAKE THE OUTPUT THE SAME AS ONE OF THE INPUTS!! */
void invmat(double *mat, double *imat)
{
    double det;
    int count;
    
    det = detmat(mat);	/* Determinant */
    adjmat(mat, imat);	/* Adjoint */
    
    for (count=0; count<9; count++)
    {
        *imat = *imat / det;
        *imat++;
    }
}


/* ====== Add two 3x3 matrices.	====== */
void addmats(double *mat1, double *mat2, double *matsum)
{
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
    *matsum++ = *mat1++ + *mat2++;
}


/* ======= Multiply two 3x3 matrices. ====== */
/*	DO NOT MAKE THE OUTPUT THE SAME AS ONE OF THE INPUTS!! */
void multmats(double *mat1, double *mat2, double *matproduct)
{
    *matproduct++ = mat1[0]*mat2[0] + mat1[3]*mat2[1] + mat1[6]*mat2[2];
    *matproduct++ = mat1[1]*mat2[0] + mat1[4]*mat2[1] + mat1[7]*mat2[2];
    *matproduct++ = mat1[2]*mat2[0] + mat1[5]*mat2[1] + mat1[8]*mat2[2];
    *matproduct++ = mat1[0]*mat2[3] + mat1[3]*mat2[4] + mat1[6]*mat2[5];
    *matproduct++ = mat1[1]*mat2[3] + mat1[4]*mat2[4] + mat1[7]*mat2[5];
    *matproduct++ = mat1[2]*mat2[3] + mat1[5]*mat2[4] + mat1[8]*mat2[5];
    *matproduct++ = mat1[0]*mat2[6] + mat1[3]*mat2[7] + mat1[6]*mat2[8];
    *matproduct++ = mat1[1]*mat2[6] + mat1[4]*mat2[7] + mat1[7]*mat2[8];
    *matproduct++ = mat1[2]*mat2[6] + mat1[5]*mat2[7] + mat1[8]*mat2[8];
}


/* Find the rotation matrix that rotates |n| radians about the vector
 * given by nx,ny,nz                                                  */
void calcrotmat(double nx, double ny, double nz, double *rmat)
{
    double ar, ai, br, bi, hp, cp, sp;
    double arar, aiai, arai2, brbr, bibi, brbi2, arbi2, aibr2, arbr2, aibi2;
    double phi;
    
    phi = sqrt(nx*nx+ny*ny+nz*nz);
    
    if (phi == 0.0) {
        *rmat++ = 1;
        *rmat++	= 0;
        *rmat++ = 0;
        *rmat++ = 0;
        *rmat++ = 1;
        *rmat++	= 0;
        *rmat++ = 0;
        *rmat++ = 0;
        *rmat++ = 1;
    } else {
        /* First define Cayley-Klein parameters */
        hp = phi/2;
        cp = cos(hp);
        sp = sin(hp)/phi;	/* /phi because n is unit length in defs. */
        ar = cp;
        ai = -nz*sp;
        br = ny*sp;
        bi = -nx*sp;
        
        /* Make auxiliary variables to speed this up */
        arar = ar*ar;
        aiai = ai*ai;
        arai2 = 2*ar*ai;
        brbr = br*br;
        bibi = bi*bi;
        brbi2 = 2*br*bi;
        arbi2 = 2*ar*bi;
        aibr2 = 2*ai*br;
        arbr2 = 2*ar*br;
        aibi2 = 2*ai*bi;
        
        /* Make rotation matrix. */
        *rmat++ = arar-aiai-brbr+bibi;
        *rmat++ = -arai2-brbi2;
        *rmat++ = -arbr2+aibi2;
        *rmat++ =  arai2-brbi2;
        *rmat++ = arar-aiai+brbr-bibi;
        *rmat++ = -aibr2-arbi2;
        *rmat++ =  arbr2+aibi2;
        *rmat++ =  arbi2-aibr2;
        *rmat++ = arar+aiai-brbr-bibi;
    }
}


/*	Set a 3x1 vector to all zeros	*/
void zerovec(double *vec)
{
    *vec++=0;
    *vec++=0;
    *vec++=0;
}


int times2intervals( double *endtimes, double *intervals, long n)
/* ------------------------------------------------------------
 * Function takes the given endtimes of intervals, and
 * returns the interval lengths in an array, assuming that
 * the first interval starts at 0.
 *
 * If the intervals are all greater than 0, then this
 * returns 1, otherwise it returns 0.
 * ------------------------------------------------------------ */
{
    int allpos;
    int count;
    double lasttime;
    
    allpos=1;
    lasttime = 0.0;
    
    for (count = 0; count < n; count++) {
        intervals[count] = endtimes[count]-lasttime;
        lasttime = endtimes[count];
        if (intervals[count] <= 0)
            allpos =0;
    }
    
    return (allpos);
}


void blochsim(double *b1real, double *b1imag,
        double *xgrad, double *ygrad, double *zgrad,
        double *tsteps, int ntime, double *e1, double *e2, double df,
        double dx, double dy, double dz,
        double *dxv, double *dyv, double *dzv,
        double *mx, double *my, double *mz,
        int mode)
/* Go through time for one df and one dx,dy,dz and one dxv,dyv,dzv. */
{
    int count;
    int tcount;
    
    double rotmat[9];
    double amat[9], bvec[3];	/* A and B propagation matrix and vector.   */
    double arot[9], brot[3];	/* A and B after rotation step.             */
    double decmat[9];           /* Decay matrix for each time step.         */
    double decvec[3];           /* Recovery vector for each time step.      */
    double rotx,roty,rotz;		/* Rotation axis coordinates.               */
    double mstart[3];
    double mfinish[3];
    double imat[9], mvec[3];
    double mcurr0[3];           /* Current magnetization before rotation.   */
    double mcurr1[3];           /* Current magnetization before decay.      */
    
    eyemat(amat);               /* A is the identity matrix.                */
    eyemat(imat);               /* I is the identity matrix.                */
    
    zerovec(bvec);
    zerovec(decvec);
    zeromat(decmat);
    
    /* Linear gradient position terms */
    double gammadx = GAMMA*dx;   /* Convert to cm•rad/Hz   */
    double gammady = GAMMA*dy;   /* Convert to cm•rad/Hz   */
    double gammadz = GAMMA*dz;   /* Convert to cm•rad/Hz   */
    
    mcurr0[0] = *mx; /* Set starting x magnetization */
    mcurr0[1] = *my; /* Set starting y magnetization */
    mcurr0[2] = *mz; /* Set starting z magnetization */
    
    for (tcount = 0; tcount < ntime; tcount++) {
        
        /*	Rotation 	*/
        // N.B. The SENSE of ROTATION was changed in code on the B Hargreaves'
        // website. Use NEW (2013) convention here, but keep z-rotation
        // following the convention in M. Levitt. "Spin Dynamics" for the
        // (observable) -1 coherence order.
        rotz = (*xgrad++ * gammadx + *ygrad++ * gammady + *zgrad++ * gammadz + df * TWOPI ) * *tsteps; // BH code is -(*xgrad++ ...)
        rotx = (+ *b1real++ * GAMMA * *tsteps ); // BH code is (- *b1real++ ...)
        roty = (+ *b1imag++ * GAMMA * *tsteps );
        // End of change.
        
        calcrotmat(rotx, roty, rotz, rotmat);
        
        if (mode == 1) {
            multmats(rotmat,amat,arot);
            multmatvec(rotmat,bvec,brot);
        } else
            multmatvec(rotmat,mcurr0,mcurr1);
        
        /* Decay */
        decvec[2]= 1- *e1;
        decmat[0]= *e2;
        decmat[4]= *e2++;
        decmat[8]= *e1++;
        
        if (mode == 1) {
            multmats(decmat,arot,amat);
            multmatvec(decmat,brot,bvec);
            addvecs(bvec,decvec,bvec);
        } else {
            multmatvec(decmat,mcurr1,mcurr0);
            addvecs(mcurr0,decvec,mcurr0);
        }
        
        /* Only do this if transient! */
        if (mode == 2) { /* Sample output at times.  */
            *mx = mcurr0[0];
            *my = mcurr0[1];
            *mz = mcurr0[2];
            mx++;
            my++;
            mz++;
        }
        
        // Update position based on velocity
        if ((*dxv!=0.0) || (*dyv!=0.0) || (*dzv!=0.0)) {
            if (*dxv!=0.0) {
                dx += *tsteps * *dxv;
                gammadx = GAMMA*dx;
            }
            if (*dyv!=0.0) {
                dy += *tsteps * *dyv;
                gammady = GAMMA*dy;
            }
            if (*dzv!=0.0) {
                dz += *tsteps++ * *dzv;
                gammadz = GAMMA*dz;
            }
        }
    } // End of time loop
    
    /* If only recording the endpoint, either store the last
       point, or calculate the steady-state endpoint. */
    if (mode==0)            /* Indicates start at given m, or m0. */
    {
        *mx = mcurr0[0];
        *my = mcurr0[1];
        *mz = mcurr0[2];
        
    } else if (mode==1)	{   /* Indicates to find steady-state magnetization */
        scalemat(amat,-1.0);		/* Negate A matrix          */
        addmats(amat,imat,amat);	/* Now amat = (I-A)         */
        invmat(amat,imat);          /* Reuse imat as inv(I-A) 	*/
        multmatvec(imat,bvec,mvec);	/* Now M = inv(I-A)*B		*/
        *mx = mvec[0];
        *my = mvec[1];
        *mz = mvec[2];
    }
}


void blochsimfz(double *b1real, double *b1imag, double *xgrad, double *ygrad, double *zgrad,
        double *tsteps, int ntime, double t1, double t2, double *dfreq, int nfreq,
        double *dxpos, double *dypos, double *dzpos, int npos,
        double *dxvel, double *dyvel, double *dzvel, int nvel,
        double *mx, double *my, double *mz, int mode)
{
    int count;
    int poscount;
    int fcount;
    int vcount;
    int totpoints;
    int totcount = 0;
    int ntout;
    
    double *e1;
    double *e2;
    double *e1ptr;
    double *e2ptr;
    double *tstepsptr;
    double *dxptr, *dyptr, *dzptr;
    double *dfptr;
    
    if (mode & 2)
        ntout = ntime;
    else
        ntout = 1;
    
    /* First calculate the E1 and E2 values at each time step. */
    e1 = (double *) malloc(ntime * sizeof(double));
    e2 = (double *) malloc(ntime * sizeof(double));
    e1ptr = e1;
    e2ptr = e2;
    tstepsptr = tsteps;
    for (count=0; count < ntime; count++) {
        *e1ptr++ = exp(- *tstepsptr / t1);
        *e2ptr++ = exp(- *tstepsptr++ / t2);
    }
    
    totpoints = npos*nfreq*nvel;
    
    for (vcount=0; vcount < nvel; vcount++) {
        dfptr = dfreq;
        
        for (fcount=0; fcount < nfreq; fcount++) {
            dxptr = dxpos;
            dyptr = dypos;
            dzptr = dzpos;
            
            for (poscount=0; poscount < npos; poscount++) {
                
                if (mode == 3) { /* Steady state AND record all time points.     */
                    /* First go through and find steady state, then */
                    /* repeat as if transient starting at steady st.*/
                    
                    blochsim(b1real, b1imag, xgrad, ygrad, zgrad,
                            tsteps, ntime, e1, e2, *dfptr, *dxptr, *dyptr, *dzptr,
                            dxvel, dyvel, dzvel, mx, my, mz, 1);
                    
                    blochsim(b1real, b1imag, xgrad, ygrad, zgrad,
                            tsteps, ntime, e1, e2, *dfptr, *dxptr++, *dyptr++, *dzptr++,
                            dxvel, dyvel, dzvel, mx, my, mz, 2);
                    
                } else {
                    blochsim(b1real, b1imag, xgrad, ygrad, zgrad,
                            tsteps, ntime, e1, e2, *dfptr, *dxptr++, *dyptr++, *dzptr++,
                            dxvel, dyvel, dzvel, mx, my, mz, mode);
                }
                
                mx += ntout;
                my += ntout;
                mz += ntout;
                
                totcount++;
                if ((totpoints > 40000) && ( ((10*totcount)/totpoints) > (10*(totcount-1)/totpoints) ))
                { DEBUG_printf("%d%% Complete.\n",(100*totcount/totpoints)); }
            }
            dfptr++;
        }
        dxvel++;
        dyvel++;
        dzvel++;
    }
    free(e1);
    free(e2);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/* bloch(b1,grad,dt,t1,t2,df,dxyz,dvxyz,mode,mx,my,mz) */
{
    double *b1r;	   /* Real-part of B1 field             */
    double *b1i;	   /* Imag-part of B1 field             */
    double *gx;        /* X              gradient (Hz/cm)   */
    double *gy;        /* Y              gradient (Hz/cm)   */
    double *gz;        /* Z              gradient (Hz/cm)   */
    double *tp;        /* Time steps (s)                    */
    double *ti;        /* Time intervals (s)                */
    double t1;         /* T1 time constant (s)              */
    double t2;         /* T2 time constant (s)              */
    double *df;        /* Off-resonance Frequencies (Hz)	*/
    double *dx;        /* X Positions (cm)                  */
    double *dy;        /* Y Positions (cm)                  */
    double *dz;        /* Z Positions (cm)                  */
    double *dxv;       /* X Velocities (cm/s)               */
    double *dyv;       /* Y Velocities (cm/s)               */
    double *dzv;       /* Z Velocities (cm/s)               */
    int md;            /* Mode - 0=from M0, 1=steady-state	*/
    double *mxin;      /* Input points                      */
    double *myin;
    double *mzin;
    double tstep;	/* Time step, if single parameter */
    double *mxout;	/* Input points  */
    double *myout;  /* Input points  */
    double *mzout;  /* Input points  */
    double *mx;	    /* Output Arrays */
    double *my;     /* Output Arrays */
    double *mz;     /* Output Arrays */
    
    int gyaflag=0;        /* 1 if gy was allocated.        */
    int gzaflag=0;        /* 1 if gz was allocated.        */
    int dyaflag=0;        /* 1 if dy was allocated.        */
    int dzaflag=0;        /* 1 if dz was allocated.        */
    
    int ntime;        /* Number of time points.              */
    int ntout;        /* Number of time poitns at output.    */
    int outsize[4];	  /* Output matrix sizes                 */
    int ngrad;        /* Number of gradient dimensions       */
    int nf;           /* Number of off-resonance frequencies */
    int npos;         /* Number of positions.  Calculated from nposN and nposM, depends on them. */
    int nposM;        /* Height of passed position matrix.   */
    int nposN;        /* Width of passed position matrix.    */
    int nfnpos;       /* Number of frequencies * number of positions. */
    int count;
    int gcount=1;     /* Gradient counter                   */
    int nvel;         /* Number of velocities.  Calculated from nvelN and nvelM. */
    int nvelM;        /* Height of passed velocity matrix.  */
    int nvelN;        /* Width of passed velocity matrix.   */
    int nfnposnvel;   /* Number of frequencies * positions * velocities.                */
    int ntnfnposnvel; /* Number of output times * frequencies * positions * velocities. */
    int dyvaflag;     /* 1 if dyv was allocated.            */
    int dzvaflag;     /* 1 if dzv was allocated.            */
    int noutdim;      /* Number of output matrix dimensions */
    
    // CTR: ERROR CHECKING. Test number of inputs.
    
    // Special case - allow bloch('debug',true) or bloch('debug',false) syntax.
    if (nrhs == 2) {
        char str1[1024];
        str1[0] = '\0';
        
        if (mxGetString(prhs[0],str1,sizeof(str1)-1)==0) {
            if (strcmp(str1,"debug") == 0) {
                debugflag = mxIsLogicalScalarTrue(prhs[1]);
                mexPrintf("Setting debug flag to %s.\n", debugflag ? "true" : "false" );
                return;
            }
        }
    }

    // Special case - allow bloch('gamma') syntax.
    if (nrhs == 1) {
        char str1[1024];
        str1[0] = '\0';
        
        if (mxGetString(prhs[0],str1,sizeof(str1)-1)==0) {
            if (strcmp(str1,"gamma") == 0) {
                if (nlhs < 1) {
                    mexPrintf("GAMMA = %g.\n", GAMMA );
                } else {
                    double *gamma_return;
                    
                    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
                    gamma_return = mxGetPr(plhs[0]);
                    *gamma_return = GAMMA;
                }
                return;
            }
        }
    }
    
    // Check number of inputs and outputs:
    if (nrhs < 8) {
        mexPrintf("Hint: Type 'doc %s' for help on input parameters.\n",mexFunctionName());
        mexErrMsgIdAndTxt("bloch:BadNInput","At least 8 inputs required.");
    }
    
    // Check all inputs are of type "double" (or char for first input):
    for (count = 0; count < nrhs; count++) {
        if (!(
                mxIsDouble(prhs[count]) // Any param can be double.
                || (count == 0 && mxIsChar(prhs[count])) // Or 1st can be char.
                ))
        {
            mexPrintf("Hint: Type 'doc %s' for help on input parameters.\n",mexFunctionName());
            mexErrMsgIdAndTxt("bloch:BadNInput","All inputs must be of type double.");
        }
    }
    // End CTR.
    
    #ifdef DEBUG
        DEBUG_printf("----------------------------------------------------------\n");
        DEBUG_printf("3D-position, 1D-frequency, and 3D-velocity Bloch Simulator\n");
        DEBUG_printf("with constant, linear, quadratic, and cubic gradients     \n");
        DEBUG_printf("----------------------------------------------------------\n\n");
    #endif
    
    ntime = mxGetM(prhs[0]) * mxGetN(prhs[0]);	/* Number of Time, RF, and Grad points */
    
    /* ====================== RF (B1) =========================
     * :  If complex, split up.  If real, allocate an imaginary part. ==== */
    if (mxIsComplex(prhs[0])) {
        b1r = mxGetPr(prhs[0]);
        b1i = mxGetPi(prhs[0]);
        
    } else {
        b1r = mxGetPr(prhs[0]);
        b1i = (double *)malloc(ntime * sizeof(double));
        for (count=0; count < ntime; count++)
            b1i[count]=0.0;
    }
    #ifdef DEBUG
        DEBUG_printf("%d B1 points.\n",ntime);
    #endif
    if (b1r == NULL)
        mexErrMsgIdAndTxt("bloch:BadPosition","B1 is not allocated.");
        
    /* ======================= Gradients ========================= */
    
    ngrad = mxGetM(prhs[1]) * mxGetN(prhs[1]);	/* Number of Time, RF, and Grad points */
    gx = mxGetPr(prhs[1]); /* X gradient is first N points. */
    #ifdef DEBUG
        DEBUG_printf("X gradient is set.");
        DEBUG_printf(" gcount = %i.\n",gcount);
    #endif
            
    if (ngrad < ++gcount*ntime) {   /* Need to allocate Y gradient. */
        #ifdef DEBUG
            DEBUG_printf("Assuming 1-dimensional gradient.\n");
        #endif
        gy = (double *)malloc(ntime * sizeof(double));
        gyaflag=1;
        for (count=0; count<ntime; count++) { gy[count]=0.0; }
    } else {
        #ifdef DEBUG
            DEBUG_printf("Y gradient is set.");
            DEBUG_printf(" gcount = %i.\n",gcount);
        #endif
        gy = gx + (gcount-1)*ntime;	/* Assign from Nx3 input array. */
    }
    
    if (ngrad < ++gcount*ntime) {  /* Need to allocate Z gradient. */
        gz = (double *)malloc(ntime * sizeof(double));
        gzaflag=1;
        for (count=0; count<ntime; count++) { gz[count]=0.0; }
    } else {
        #ifdef DEBUG
            DEBUG_printf("Z gradient is set.");
            DEBUG_printf(" gcount = %i.\n",gcount);
        #endif
        gz = gx + (gcount-1)*ntime; /* Assign from Nx3 input array. */
    }
    
    /* Warning if Gradient length is not an integer multiple of the RF length. */
    #ifdef DEBUG
        DEBUG_printf("%d Gradient Points (total).\n",ngrad);
    #endif
    /* if ( (ngrad != ntime) && (ngrad != 2*ntime) && (ngrad != 3*ntime) ) */
    if ( (ngrad % ntime) > 0 )
        mexErrMsgIdAndTxt("bloch:BadGradientLength","Gradient length differs from B1 length.");
    if (gx == NULL)
        mexErrMsgIdAndTxt("bloch:BadGradientLength","gx is not allocated.");
    if (gy == NULL)
        mexErrMsgIdAndTxt("bloch:BadGradientLength","gy is not allocated.");
    if (gz == NULL)
        mexErrMsgIdAndTxt("bloch:BadGradientLength","gz is not allocated.");
            
    /* === Time points ===== */
        
    /*	THREE Cases:
        1) Single value given -> this is the interval length for all.
        2) List of intervals given.
        3) Monotonically INCREASING list of end times given.
	For all cases, the goal is for tp to have the intervals. */
        
    ti = NULL;
    tp = mxGetPr(prhs[2]);
    if (tp == NULL)
        mexErrMsgIdAndTxt("bloch:BadPosition","tp is not allocated.");
    
    if (mxGetM(prhs[2]) * mxGetN(prhs[2]) == 1)	{ /* === Case 1 === */
        tp = (double *)malloc(ntime * sizeof(double));
        tstep = *(mxGetPr(prhs[2]));
        for (count =0; count < ntime; count++)
            tp[count]=tstep;
        
    } else if (mxGetM(prhs[2]) * mxGetN(prhs[2]) != ntime)
        mexErrMsgIdAndTxt("bloch:BadB1Length","Time-point length differs from B1 length.");
    
    else {
        tp = mxGetPr(prhs[2]);
        ti = (double *)malloc(ntime * sizeof(double));
        if (( times2intervals( tp, ti, ntime ))) {
            DEBUG_printf("Times are monotonically increasing.\n");
            tp = ti;
        }
    }
    
    /* === Relaxation Times ===== */
    if (mxGetPr(prhs[3]) == NULL)
        mexErrMsgIdAndTxt("bloch:BadPosition","t1 is not allocated.");
    if (mxGetPr(prhs[4]) == NULL)
        mexErrMsgIdAndTxt("bloch:BadPosition","t2 is not allocated.");
    
    t1 = *mxGetPr(prhs[3]);
    t2 = *mxGetPr(prhs[4]);
    
    #ifdef DEBUG
        DEBUG_printf("t1 = %d \n",t1);
        DEBUG_printf("t2 = %d \n",t2);
    #endif
        
    /* === Frequency Points ===== */
    df = mxGetPr(prhs[5]);
    nf = mxGetM(prhs[5]) * mxGetN(prhs[5]);
    
    #ifdef DEBUG
        DEBUG_printf("%d Frequency points.\n",nf);
    #endif
    if (df == NULL)
        mexErrMsgIdAndTxt("bloch:BadPosition","df is not allocated.");

    /* === Position Points ===== */
    nposM = mxGetM(prhs[6]);
    nposN = mxGetN(prhs[6]);
    
    #ifdef DEBUG
        DEBUG_printf("Position vector is %d x %d.\n",nposM,nposN);
    #endif
    
    if (nposN==3) { /* Assume 3 position dimensions given */
        npos = nposM;
        #ifdef DEBUG
            DEBUG_printf("Assuming %d 3-Dimensional Positions.\n",npos);
        #endif
        dx = mxGetPr(prhs[6]);
        dy = dx + npos;
        dz = dy + npos;

    } else if (nposN==2) { /* Assume only 2 position dimensions given */
        npos = nposM;
        #ifdef DEBUG
            DEBUG_printf("Assuming %d 2-Dimensional Positions.\n",npos);
        #endif
        dx = mxGetPr(prhs[6]);
        dy = dx + npos;
        dz = (double *)malloc(npos * sizeof(double));
        dzaflag=1;
        for (count=0; count < npos; count++)
            dz[count]=0.0;
    
    } else { /* Either 1xN, Nx1 or something random.  In all these
     * cases we assume that 1 position is given, because it
     * is too much work to try to figure out anything else! */
        npos = nposM * nposN;
        #ifdef DEBUG
            DEBUG_printf("Assuming %d 1-Dimensional Positions.\n",npos);
        #endif
        dx = mxGetPr(prhs[6]);
        dy = (double *)malloc(npos * sizeof(double));
        dz = (double *)malloc(npos * sizeof(double));
        dyaflag=1;
        dzaflag=1;
        for (count=0; count < npos; count++) {
            dy[count]=0.0;
            dz[count]=0.0;
        }
        #ifdef DEBUG
            if ((nposM !=1) && (nposN!=1)) {
                DEBUG_printf("Position vector should be 1xN, Nx1, Nx2 or Nx3.\n");
                DEBUG_printf(" -> Assuming 1 position dimension is given.\n");
            }
        #endif
    }
    if (dx == NULL)
        mexErrMsgIdAndTxt("bloch:BadPosition","dx is not allocated.");
    if (dy == NULL)
        mexErrMsgIdAndTxt("bloch:BadPosition","dy is not allocated.");
    if (dz == NULL)
        mexErrMsgIdAndTxt("bloch:BadPosition","dz is not allocated.");

    nfnpos = nf*npos;	/* Just used to speed things up below. 	*/ 

    /* === Velocity Points ===== */
    nvelM = mxGetM(prhs[7]);
    nvelN = mxGetN(prhs[7]);
    #ifdef DEBUG
        DEBUG_printf("Velocity vector is %d x %d.\n",nvelM,nvelN);
	#endif

    if (nvelN==3) {         /* Assume 3 velocity dimensions given */
        nvel = nvelM;
        #ifdef DEBUG
            DEBUG_printf("Assuming %d 3-Dimensional Velocities.\n",nvel);
        #endif
        dxv = mxGetPr(prhs[7]);
        dyv = dxv + nvel;
        dzv = dyv + nvel;
    } else if (nvelN==2) {	/* Assume only 2 velocity dimensions given */
        nvel = nvelM;
        #ifdef DEBUG
            DEBUG_printf("Assuming %d 2-Dimensional Velocities.\n",nvel);
        #endif
        dxv = mxGetPr(prhs[7]);
        dyv = dxv + nvel;
        dzv = (double *)malloc(nvel * sizeof(double));
        dzvaflag=1;
        for (count=0; count < nvel; count++)
            dzv[count]=0.0;
    } else {               	/* Assume only 1 velocity dimension given */
        nvel = nvelM * nvelN;
        #ifdef DEBUG
            DEBUG_printf("Assuming %d 1-Dimensional Velocities.\n",npos);
        #endif
        dxv = mxGetPr(prhs[7]);
        dyv = (double *)malloc(nvel * sizeof(double));
        dzv = (double *)malloc(nvel * sizeof(double));
        dyvaflag=1;
        dzvaflag=1;
        for (count=0; count < nvel; count++) {
            dyv[count]=0.0;
            dzv[count]=0.0;
        }
        #ifdef DEBUG
            if ((nvelM !=1) && (nvelN!=1)) {
                DEBUG_printf("Velocity vector should be 1xN, Nx1, Nx2 or Nx3.\n");
                DEBUG_printf(" -> Assuming 1 velocity dimension is given.\n");
            }
        #endif
    }
    if (dxv == NULL)
        mexErrMsgIdAndTxt("bloch:BadVelocity","dxv is not allocated.");
    if (dyv == NULL)
        mexErrMsgIdAndTxt("bloch:BadVelocity","dyv is not allocated.");
    if (dzv == NULL)
        mexErrMsgIdAndTxt("bloch:BadVelocity","dzv is not allocated.");

    nfnposnvel = nfnpos * nvel;
  
    /* ===== Mode, defaults to 0 (simulate single endpoint, transient). ==== */
    if (nrhs > 8)
        md = (int)(*mxGetPr(prhs[8]));
    else
        md = 0;

    if (md & 2)
        ntout = ntime;		/* Include time points.	*/
    else
        ntout = 1;

    #ifdef DEBUG
        DEBUG_printf("Mode = %d, %d Output Time Points.\n",md,ntout);
    #endif

    ntnfnposnvel = ntout*nfnposnvel;

    #ifdef DEBUG
        if ((md & 1)==0)
            DEBUG_printf("Simulation from Initial Condition.\n");
        else
            DEBUG_printf("Simulation of Steady-State.\n");
    
        if ((md & 2)==0)
            DEBUG_printf("Simulation to Endpoint.\n");
        else
            DEBUG_printf("Simulation over Time.\n");
    #endif

    /* ===== Allocate Output Magnetization vectors arrays.	*/
    plhs[0] = mxCreateDoubleMatrix(ntnfnposnvel,1,mxREAL);	/* Mx, output. */
    plhs[1] = mxCreateDoubleMatrix(ntnfnposnvel,1,mxREAL);	/* My, output. */
    plhs[2] = mxCreateDoubleMatrix(ntnfnposnvel,1,mxREAL);	/* Mz, output. */
    
    mx = mxGetPr(plhs[0]);
    my = mxGetPr(plhs[1]);
    mz = mxGetPr(plhs[2]);
    
    mxout = mx;
    myout = my;
    mzout = mz;

    /* ===== If Initial Magnetization is given... */
    if ( (nrhs > 11) &&
            (mxGetM(prhs[9])  * mxGetN(prhs[9])  == nfnposnvel) &&
            (mxGetM(prhs[10]) * mxGetN(prhs[10]) == nfnposnvel) &&
            (mxGetM(prhs[11]) * mxGetN(prhs[11]) == nfnposnvel)  )
        
        /* Set output magnetization to that passed.
         * If multiple time points, then just the
         * first is set.                            */
    {
    #ifdef DEBUG
        DEBUG_printf("Using Specified Initial Magnetization.\n");
    #endif
    
    mxin = mxGetPr(prhs[9]);
    myin = mxGetPr(prhs[10]);
    mzin = mxGetPr(prhs[11]);
    for (count =0; count < nfnposnvel; count++)
    {
        *mxout = *mxin++;
        *myout = *myin++;
        *mzout = *mzin++;
        mxout += ntout;
        myout += ntout;
        mzout += ntout;
    }
    } else {
        #ifdef DEBUG
            if (nrhs > 11) { /* Magnetization given, but wrong size! */
                mexErrMsgIdAndTxt("bloch:BadMagnetization","Initial magnetization passed, but not Npositions x Nfreq x Nvelocities.");
            }
            DEBUG_printf(" --> Using [0; 0; 1] for initial magnetization.\n");
        #endif
        for (count =0; count < nfnposnvel; count++) {
            *mxout = 0;	/* Set magnetization to Equilibrium */
            *myout = 0;
            *mzout = 1;
            mxout += ntout;
            myout += ntout;
            mzout += ntout;
        }
    }


    /* ======= Do The Simulation! ====== */
    #ifdef DEBUG
        DEBUG_printf("Calling blochsimfz() function in Mex file.\n");
    #endif

    blochsimfz(b1r,b1i,gx,gy,gz,tp,ntime,t1,t2,df,nf,dx,dy,dz,npos,dxv,dyv,dzv,nvel,mx,my,mz,md);

    
    /* ======= Reshape Output Matrices ====== */
    noutdim = (int)(ntout>1) + (int)(npos>1) + (int)(nf>1) + (int)(nvel>1);

    if (noutdim == 4) {
        outsize[0]=ntout;
        outsize[1]=npos;
        outsize[2]=nf;
        outsize[3]=nvel;
        mxSetDimensions(plhs[0],outsize,4);  /* Set to 4D array. */
        mxSetDimensions(plhs[1],outsize,4);  /* Set to 4D array. */
        mxSetDimensions(plhs[2],outsize,4);  /* Set to 4D array. */
        mxSetDimensions(plhs[3],outsize,4);  /* Set to 4D array. */
    } else if (noutdim == 3) { /* Try 3 dimensions */
        if ((ntout > 1) && (npos > 1) && (nf > 1) && (nvel == 1)) {
            outsize[0]=ntout;
            outsize[1]=npos;
            outsize[2]=nf;
        } else if ((ntout > 1) && (npos > 1) && (nf == 1) && (nvel > 1)) {
            outsize[0]=ntout;
            outsize[1]=npos;
            outsize[2]=nvel;
        } else if ((ntout > 1) && (npos == 1) && (nf > 1) && (nvel > 1)) {
            outsize[0]=ntout;
            outsize[1]=nf;
            outsize[2]=nvel;
        } else if ((ntout == 1) && (npos > 1) && (nf > 1) && (nvel > 1)) {
            outsize[0]=npos;
            outsize[1]=nf;
            outsize[2]=nvel;
        }
        mxSetDimensions(plhs[0],outsize,3);  /* Set to 3D array. */
        mxSetDimensions(plhs[1],outsize,3);  /* Set to 3D array. */
        mxSetDimensions(plhs[2],outsize,3);  /* Set to 3D array. */
    } else { /* Only 2 or 1 dimensions */
        if (ntout > 1) {
            outsize[0]=ntout;
            outsize[1]=npos*nf*nvel;
        } else if (npos > 1) {
            outsize[0]=npos;
            outsize[1]=nf*nvel;
        } else {
            outsize[0]=nf;
            outsize[1]=nvel;
        }
        mxSetDimensions(plhs[0],outsize,2);  /* Set to 2D array. */
        mxSetDimensions(plhs[1],outsize,2);  /* Set to 2D array. */
        mxSetDimensions(plhs[2],outsize,2);  /* Set to 2D array. */
    }

    /* ====== Free up allocated memory, if necessary. ===== */
    if (!mxIsComplex(prhs[0])) free(b1i);
    if (mxGetM(prhs[2]) * mxGetN(prhs[2]) == 1) free(tp);
    if (ti != NULL) free(ti);
    if (dyaflag==1) free(dy);
    if (dzaflag==1) free(dz);
    if (gyaflag==1) free(gy);
    if (gzaflag==1) free(gz);
    if (dyvaflag==1) free(dyv);
    if (dzvaflag==1) free(dzv);
}
